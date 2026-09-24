"""HTML/PDF reports: motifs.html table handling (MEME tomtom q-values and TOMTOM-lite p-values) and report mains."""
import argparse
import json
import os

import matplotlib

matplotlib.use("Agg")

import pandas as pd  # noqa: E402
import pytest  # noqa: E402
from modiscolite.report import path_to_image_html  # noqa: E402

from chrombpnet.evaluation.modisco import convert_html_to_pdf  # noqa: E402
from chrombpnet.helpers.generate_reports import make_html, make_html_bias, modisco_table  # noqa: E402

PATTERNS = [("pos_patterns", 0, 812), ("pos_patterns", 1, 344), ("pos_patterns", 2, 97),
            ("neg_patterns", 0, 250), ("neg_patterns", 1, 41)]


def write_motifs_html(path, stat="qval", n_matches=3, patterns=PATTERNS, names=("CTCF_HUMAN", "GATA1_MA0035.4",
                                                                               "Tn5_bias_1")):
    """motifs.html exactly as modiscolite.report.report_motifs writes it (same columns, formatters, flags)."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    rows = []
    for k, (group, idx, n) in enumerate(patterns):
        tag = "{}.pattern_{}".format(group, idx)
        row = {"pattern": tag, "num_seqlets": n,
               "modisco_cwm_fwd": "./trimmed_logos/{}.cwm.fwd.png".format(tag),
               "modisco_cwm_rev": "./trimmed_logos/{}.cwm.rev.png".format(tag)}
        for i in range(n_matches):
            match = None if (stat == "qval" and k == len(patterns) - 1 and i > 0) else names[(k + i) % len(names)]
            row["match{}".format(i)] = match
            row["{}{}".format(stat, i)] = None if match is None else 10.0 ** -(k + i + 3)
            row["match{}_logo".format(i)] = "NA" if match is None else "./{}.png".format(match)
        rows.append(row)
    columns = ["pattern", "num_seqlets", "modisco_cwm_fwd", "modisco_cwm_rev"]
    for i in range(n_matches):
        columns += ["match{}".format(i), "{}{}".format(stat, i), "match{}_logo".format(i)]
    df = pd.DataFrame(rows)[columns]
    formatters = dict(modisco_cwm_fwd=path_to_image_html, modisco_cwm_rev=path_to_image_html)
    for i in range(n_matches):
        formatters["match{}_logo".format(i)] = path_to_image_html
    with open(path, "w") as f:
        df.to_html(f, escape=False, formatters=formatters, index=False)
    return path


def write_plain_motifs_html(path):
    """motifs.html written without a motif database (no TOMTOM columns)."""
    return write_motifs_html(path, n_matches=0)


@pytest.mark.parametrize("stat", ["qval", "pval"])
def test_match_stat_and_columns(tmp_path, stat):
    html = open(write_motifs_html(str(tmp_path / "motifs.html"), stat=stat)).read()
    assert modisco_table.match_stat(html) == stat
    assert modisco_table.table_columns(html) == [
        "pattern", "num_seqlets", "modisco_cwm_fwd", "modisco_cwm_rev",
        "match0", stat + "0", "match0_logo", "match1", stat + "1", "match1_logo", "match2", stat + "2", "match2_logo"]


def test_match_stat_without_motif_database(tmp_path):
    html = open(write_plain_motifs_html(str(tmp_path / "motifs.html"))).read()
    assert modisco_table.match_stat(html) is None
    assert modisco_table.stat_description(None) == "qvals (qval0,qval1,qval2)"


def test_stat_wording():
    assert modisco_table.stat_description("qval") == "qvals (qval0,qval1,qval2)"
    assert modisco_table.stat_description("pval").startswith("p-values (pval0,pval1,pval2; TOMTOM-lite")
    assert modisco_table.stat_name("qval") == "qvals" and modisco_table.stat_name("pval") == "p-values"


@pytest.mark.parametrize("stat,n_matches", [("qval", 3), ("pval", 3), ("qval", 5), ("pval", 1), ("qval", 0)])
def test_drop_negative_patterns_for_any_column_layout(tmp_path, stat, n_matches):
    # the chrombpnet 1.x parser assumed 15 html lines per row (exactly 3 matches); this must not
    html = open(write_motifs_html(str(tmp_path / "motifs.html"), stat=stat, n_matches=n_matches)).read()
    kept = modisco_table.drop_negative_patterns(html)
    assert "neg_patterns" not in kept
    for tag in ("pos_patterns.pattern_0", "pos_patterns.pattern_1", "pos_patterns.pattern_2"):
        assert kept.count("<td>{}</td>".format(tag)) == 1
    assert kept.count("<tr") == 1 + 3  # header + positive patterns
    assert kept.count("<td") == 3 * (4 + 3 * n_matches)
    assert kept.rstrip().endswith("</table>")
    # the remaining rows are untouched, byte for byte
    body = html[html.index("<tbody>"):html.index("<td>neg_patterns")]
    assert body[:body.rindex("</tr>")] in kept


def test_drop_negative_patterns_without_negatives(tmp_path):
    html = open(write_motifs_html(str(tmp_path / "motifs.html"), patterns=PATTERNS[:2])).read()
    assert modisco_table.drop_negative_patterns(html) == html


@pytest.mark.parametrize("stat", ["qval", "pval"])
def test_load_motifs_table(tmp_path, stat):
    path = write_motifs_html(str(tmp_path / "modisco_counts" / "motifs.html"), stat=stat)
    table, kind = modisco_table.load_motifs_table(path, "modisco_counts", drop_negative=True)
    assert kind == stat
    assert table.startswith('<table border="1" class="new">')
    assert '<img src="./modisco_counts/trimmed_logos/pos_patterns.pattern_1.cwm.fwd.png" width="240", class="cover" >' \
        in table
    assert '<img src="./modisco_counts/CTCF_HUMAN.png" width="240", class="cover" >' in table
    assert ">pos__0<" in table and ">pos__2<" in table and "neg_" not in table
    assert "<th>NumSeqs</th>" in table and "<th>cwm_fwd</th>" in table and "<th>cwm_rev</th>" in table
    assert "<th>{}0</th>".format(stat) in table

    table, _ = modisco_table.load_motifs_table(path, "modisco_counts")
    assert ">neg__0<" in table and ">neg__1<" in table


# ---------------------------------------------------------------- report mains

def weasyprint_or_skip():
    try:
        import weasyprint  # noqa: F401
    except (ImportError, OSError) as e:
        pytest.skip("weasyprint unavailable: {}".format(e))


def metrics(pearsonr):
    split = {"pearsonr": pearsonr, "spearmanr": pearsonr - 0.05, "mse": 1.25}
    prof = {"median_jsd": 0.61, "median_norm_jsd": 0.18}
    return {"counts_metrics": {"peaks": split, "nonpeaks": split, "peaks_and_nonpeaks": split},
            "profile_metrics": {"peaks": prof, "nonpeaks": prof, "peaks_and_nonpeaks": prof}}


def write_training_log(path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    pd.DataFrame({"epoch": [0, 1, 2], "loss": [310.2, 300.1, 298.7], "val_loss": [305.0, 301.4, 300.9],
                  "logcount_predictions_loss": [2.1, 1.9, 1.8],
                  "logits_profile_predictions_loss": [300.0, 290.0, 288.0]}).to_csv(path, index=False)


def is_pdf(path):
    with open(path, "rb") as f:
        return f.read(5) == b"%PDF-"


@pytest.mark.parametrize("stat", ["qval", "pval"])
def test_make_html_bias_pipeline_report(tmp_path, stat):
    weasyprint_or_skip()
    out = tmp_path / "bias_model"
    write_training_log(str(out / "logs" / "bias.log"))
    (out / "evaluation").mkdir(parents=True, exist_ok=True)
    json.dump(metrics(0.12), open(str(out / "evaluation" / "bias_metrics.json"), "w"))
    write_motifs_html(str(out / "evaluation" / "modisco_profile" / "motifs.html"), stat=stat)
    write_motifs_html(str(out / "evaluation" / "modisco_counts" / "motifs.html"), stat=stat, n_matches=2)

    args = argparse.Namespace(input_dir=str(out), file_prefix=None, command="pipeline", html_prefix="./")
    make_html_bias.main(args)

    html = (out / "evaluation" / "overall_report.html").read_text()
    assert is_pdf(str(out / "evaluation" / "overall_report.pdf"))
    assert (out / "evaluation" / "epoch_loss.png").stat().st_size > 0
    assert "neg_" not in html and ">pos__0<" in html
    assert "./modisco_profile/trimmed_logos/" in html and "./modisco_counts/trimmed_logos/" in html
    if stat == "qval":
        assert "The qvals (qval0,qval1,qval2) should be high" in html
        assert "The qvals should be low if the closest hit is enzyme bias motif" in html
        assert "p-values" not in html
    else:
        assert "The p-values (pval0,pval1,pval2; TOMTOM-lite" in html
        assert "The p-values should be low if the closest hit" in html
        assert "qval" not in html


def test_make_html_bias_train_report_with_prefix(tmp_path):
    weasyprint_or_skip()
    out = tmp_path / "bias_model"
    write_training_log(str(out / "logs" / "sample_bias.log"))
    (out / "evaluation").mkdir(parents=True, exist_ok=True)
    args = argparse.Namespace(input_dir=str(out), file_prefix="sample", command="train", html_prefix="https://x/")
    make_html_bias.main(args)
    html = (out / "evaluation" / "sample_overall_report.html").read_text()
    assert "https://x/sample_bw_shift_qc.png" in html  # html_prefix applied to the final html
    assert is_pdf(str(out / "evaluation" / "sample_overall_report.pdf"))


@pytest.mark.parametrize("stat,data_type", [("qval", "ATAC"), ("pval", "DNASE")])
def test_make_html_pipeline_report(tmp_path, stat, data_type):
    weasyprint_or_skip()
    out = tmp_path / "chrombpnet_model"
    write_training_log(str(out / "logs" / "chrombpnet.log"))
    (out / "evaluation").mkdir(parents=True, exist_ok=True)
    json.dump(metrics(-0.1), open(str(out / "evaluation" / "bias_metrics.json"), "w"))
    json.dump(metrics(0.7), open(str(out / "evaluation" / "chrombpnet_metrics.json"), "w"))
    (out / "evaluation" / "chrombpnet_nobias_max_bias_response.txt").write_text("corrected_0.0012")
    # a brace in a motif name must not break the str.format templating of the report
    write_motifs_html(str(out / "evaluation" / "modisco_profile" / "motifs.html"), stat=stat,
                      names=("CTCF_HUMAN", "ODD{0}_MOTIF", "GATA1"))

    args = argparse.Namespace(input_dir=str(out), data_type=data_type, file_prefix=None, command="pipeline",
                              html_prefix="./")
    make_html.main(args)

    html = (out / "evaluation" / "overall_report.html").read_text()
    assert is_pdf(str(out / "evaluation" / "overall_report.pdf"))
    assert "ODD{0}_MOTIF" in html
    assert ">neg__0<" in html  # the chrombpnet report keeps negative patterns
    assert ("tn5_1.footprint.png" in html) == (data_type == "ATAC")
    if stat == "qval":
        assert "The qvals (qval0,qval1,qval2) should be low (< 0.0001)" in html
    else:
        assert "The p-values (pval0,pval1,pval2; TOMTOM-lite, not corrected for multiple testing) should be low" in html


def test_convert_html_to_pdf(tmp_path):
    weasyprint_or_skip()
    html = tmp_path / "motifs.html"
    write_motifs_html(str(html), stat="pval")
    pdf = tmp_path / "motifs.pdf"
    convert_html_to_pdf.main(str(html), str(pdf))
    assert is_pdf(str(pdf))
