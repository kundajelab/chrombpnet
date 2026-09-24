"""MoDISco / report parity against real chrombpnet 1.x outputs (modisco-lite 2.0.7 era).

Set CHROMBPNET_LEGACY_RESULTS to one or more (os.pathsep-separated) chrombpnet 1.x results or run directories;
they are searched recursively for
  - evaluation/modisco_{profile,counts}/motifs.html: the report table handling must keep exactly the rows the
    chrombpnet 1.x parser kept (fast);
  - auxiliary/interpret_subsample/<prefix>{bias,chrombpnet_nobias}.profile_scores.h5 with the sibling
    <prefix>modisco_results_profile_scores.h5: `modisco motifs` 2.5 with chrombpnet's explicit settings is rerun
    on the same scores and compared with the legacy patterns (slow: hours of CPU at 50,000 seqlets; set
    CHROMBPNET_MODISCO_THREADS to the CPUs you may use, CHROMBPNET_MODISCO_PARITY_MAX (default 1) to cap the runs).
"""
import glob
import os
import re

import h5py
import numpy as np
import pytest

from chrombpnet.evaluation.modisco import run
from chrombpnet.helpers.generate_reports import modisco_table


def roots():
    value = os.environ.get("CHROMBPNET_LEGACY_RESULTS")
    if not value:
        pytest.skip("CHROMBPNET_LEGACY_RESULTS not set")
    paths = [p for p in value.split(os.pathsep) if os.path.isdir(p)]
    if not paths:
        pytest.skip("CHROMBPNET_LEGACY_RESULTS has no existing directory: " + value)
    return paths


def find(pattern):
    return sorted({f for root in roots() for f in glob.glob(os.path.join(root, "**", pattern), recursive=True)})


def legacy_remove_negs(tables):
	# chrombpnet 1.x make_html_bias.remove_negs, verbatim
	new_lines=[]
	set_flag = True
	lines = tables.split("\n")
	jdx=0
	for idx in range(len(lines)-1):

		if jdx==15:
			set_flag = True
			jdx=0

		if "neg_" in lines[idx+1]:
			set_flag=False
			jdx = 0

		if set_flag:
			new_lines.append(lines[idx])
		else:
			jdx+=1
	new_lines.append(lines[-1])
	return  "\n".join(new_lines)


def row_names(html):
    return re.findall(r"<td>\s*((?:pos|neg)_\S*?)\s*</td>", html)


def test_negative_rows_match_legacy_parser_on_legacy_reports():
    reports = find(os.path.join("evaluation", "modisco_*", "motifs.html"))
    if not reports:
        pytest.skip("no evaluation/modisco_*/motifs.html under CHROMBPNET_LEGACY_RESULTS")
    for path in reports:
        html = open(path).read()
        assert modisco_table.match_stat(html) in ("qval", None), path
        # legacy flow: rename first, then drop rows by counting lines
        legacy_input = html.replace(">pos_patterns.pattern", ">pos_").replace(">neg_patterns.pattern", ">neg_")
        legacy = legacy_remove_negs(legacy_input)
        ours, _ = modisco_table.load_motifs_table(path, "modisco_x", drop_negative=True)
        assert row_names(ours) == row_names(legacy), path
        assert not any(name.startswith("neg_") for name in row_names(ours)), path
        assert ours.count("<img") == legacy.count("<img"), path


def score_result_pairs():
    pairs = []
    for scores in find(os.path.join("auxiliary", "interpret_subsample", "*.profile_scores.h5")):
        name = os.path.basename(scores)
        prefix = re.sub(r"(bias|chrombpnet_nobias)\.profile_scores\.h5$", "", name)
        legacy = os.path.join(os.path.dirname(scores), prefix + "modisco_results_profile_scores.h5")
        if prefix != name and os.path.isfile(legacy):
            pairs.append((scores, legacy))
    return pairs


def patterns(path):
    out = {}
    with h5py.File(path, "r") as f:
        for group in ("pos_patterns", "neg_patterns"):
            if group not in f:
                continue
            for name in f[group]:
                p = f[group][name]
                out[(group, int(name.split("_")[-1]))] = (np.array(p["contrib_scores"][:]),
                                                          int(p["seqlets"]["n_seqlets"][:][0]))
    return out


def cwm_similarity(a, b, max_offset=5):
    """Best Pearson correlation of two CWMs over small offsets and both strands."""
    best = -1.0
    for bb in (b, b[::-1, ::-1]):
        for off in range(-max_offset, max_offset + 1):
            x = a[max(0, off):len(a) + min(0, off)]
            y = bb[max(0, -off):len(bb) + min(0, -off)]
            n = min(len(x), len(y))
            if n < 10:
                continue
            best = max(best, float(np.corrcoef(x[:n].ravel(), y[:n].ravel())[0, 1]))
    return best


@pytest.mark.slow
def test_modisco_25_reproduces_legacy_patterns(tmp_path):
    pairs = score_result_pairs()
    if not pairs:
        pytest.skip("no interpret_subsample/*.profile_scores.h5 with legacy modisco results under CHROMBPNET_LEGACY_RESULTS")
    threads = os.environ.get("CHROMBPNET_MODISCO_THREADS")
    for i, (scores, legacy_h5) in enumerate(pairs[:int(os.environ.get("CHROMBPNET_MODISCO_PARITY_MAX", "1"))]):
        new_h5 = run.modisco_motifs(scores, str(tmp_path / "modisco_{}.h5".format(i)),
                                    threads=int(threads) if threads else None)
        legacy, new = patterns(legacy_h5), patterns(new_h5)
        for group in ("pos_patterns", "neg_patterns"):
            old_g = {k: v for k, v in legacy.items() if k[0] == group}
            new_g = {k: v for k, v in new.items() if k[0] == group}
            assert abs(len(new_g) - len(old_g)) <= max(3, 0.3 * len(old_g)), (scores, group, len(old_g), len(new_g))
            widths_old = {v[0].shape[0] for v in old_g.values()}
            widths_new = {v[0].shape[0] for v in new_g.values()}
            assert widths_new <= widths_old | {30}, (scores, group, widths_old, widths_new)
            top = sorted(old_g.values(), key=lambda v: -v[1])[:10]
            matched = sum(max((cwm_similarity(o[0], n[0]) for n in new_g.values()), default=-1) >= 0.9 for o in top)
            assert matched >= 0.8 * len(top), (scores, group, matched, len(top))
