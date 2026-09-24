"""Read the motifs.html table written by `modisco report-simple` for the ChromBPNet HTML/PDF reports.

The table is pandas `to_html` output with one <tr> per pattern. Its <img> cells (CWM and TOMTOM match
logos) are the content of the report, so it is edited row-wise as HTML rather than parsed to text (which
is what pandas.read_html would do). The TOMTOM statistic columns are qval0..N with MEME `tomtom` and
pval0..N with `--tomtom-lite`.
"""
import re

_HEADER_CELL = re.compile(r"<th(?:\s[^>]*)?>(.*?)</th>", re.S)
_ROW = re.compile(r"<tr(?:\s[^>]*)?>(.*?)</tr>\s*", re.S)
_FIRST_CELL = re.compile(r"\s*<td(?:\s[^>]*)?>\s*(.*?)\s*</td>", re.S)


def table_columns(html):
	"""Column names of the table header, in order."""
	return [re.sub(r"<[^>]+>", "", c).strip() for c in _HEADER_CELL.findall(html)]


def match_stat(html):
	"""'qval' (MEME tomtom), 'pval' (TOMTOM-lite) or None (no motif database given to the report)."""
	columns = table_columns(html)
	for kind in ("qval", "pval"):
		if any(re.fullmatch(kind + r"\d+", c) for c in columns):
			return kind
	return None


def drop_negative_patterns(html):
	"""Remove the rows of neg_patterns.* whatever the number of columns."""
	def keep(row):
		first = _FIRST_CELL.match(row.group(1))
		if first is not None and first.group(1).startswith("neg_"):
			return ""
		return row.group(0)
	return _ROW.sub(keep, html)


def stat_name(kind):
	return "p-values" if kind == "pval" else "qvals"


def stat_description(kind, n=3):
	"""Wording for the TOMTOM statistic columns (match_0..match_{n-1}) of a report table."""
	if kind == "pval":
		return "p-values ({}; TOMTOM-lite, not corrected for multiple testing)".format(
			",".join("pval{}".format(i) for i in range(n)))
	return "qvals ({})".format(",".join("qval{}".format(i) for i in range(n)))


def load_motifs_table(motifs_html, subdir, drop_negative=False):
	"""motifs.html as a report table (image paths rebased to ./<subdir>/, short pattern and column names)
	and the kind of its TOMTOM statistic columns (see match_stat)."""
	with open(motifs_html) as f:
		html = f.read()
	kind = match_stat(html)
	if drop_negative:
		html = drop_negative_patterns(html)
	html = html.replace("./","./{}/".format(subdir)).replace("width=\"240\"","width=\"240\", class=\"cover\"").replace(">pos_patterns.pattern",">pos_").replace(">neg_patterns.pattern",">neg_").replace("modisco_cwm_fwd","cwm_fwd").replace("modisco_cwm_rev","cwm_rev").replace("num_seqlets","NumSeqs").replace("dataframe","new")
	return html, kind
