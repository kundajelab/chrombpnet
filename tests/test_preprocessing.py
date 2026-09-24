"""Preprocessing helpers: seeded read sampling, shift decisions, subprocess error propagation, gzipped inputs,
bigWig building."""
import argparse
import gzip
import shutil
import subprocess
import sys
import warnings

import numpy as np
import pyBigWig
import pytest

from chrombpnet.data import DefaultDataFile, get_default_data_path
from chrombpnet.helpers.preprocessing import auto_shift_detect, reads_to_bigwig

CHROM_LEN = 5000


@pytest.fixture
def genome(tmp_path):
    rng = np.random.RandomState(0)
    fasta = tmp_path / "genome.fa"
    sizes = tmp_path / "genome.chrom.sizes"
    with open(fasta, "w") as f:
        for chrom in ("chr1", "chr2"):
            seq = "".join(rng.choice(list("ACGT"), CHROM_LEN))
            f.write(">{}\n".format(chrom))
            for i in range(0, CHROM_LEN, 60):
                f.write(seq[i:i + 60] + "\n")
    sizes.write_text("chr1\t{0}\nchr2\t{0}\n".format(CHROM_LEN))
    return str(fasta), str(sizes)


def write_tagalign(path, reads):
    with open(path, "w") as f:
        for chrom, start, end, strand in reads:
            f.write("{}\t{}\t{}\tN\t1000\t{}\n".format(chrom, start, end, strand))
    return str(path)


def random_reads(n, seed=1):
    rng = np.random.RandomState(seed)
    reads = []
    for i in range(n):
        chrom = "chr1" if i % 2 else "chr2"
        start = int(rng.randint(100, CHROM_LEN - 200))
        reads.append((chrom, start, start + 50, "+" if rng.rand() < 0.5 else "-"))
    return reads


def write_fragments(path, fragments):
    # chr, start, end, barcode, count (cellranger / ArchR layout)
    with open(path, "w") as f:
        for i, (chrom, start, end) in enumerate(fragments):
            f.write("{}\t{}\t{}\tBC{}\t1\n".format(chrom, start, end, i % 7))
    return str(path)


def random_fragments(n, seed=2):
    return [(chrom, start, end) for chrom, start, end, _ in random_reads(n, seed)]


def write_input(tmp_path, kind, n, seed=3):
    """A plain fragment or tagAlign file and the reads (chr, start, end, strand) it stands for."""
    if kind == "fragment":
        fragments = random_fragments(n, seed)
        path = write_fragments(tmp_path / "reads.tsv", fragments)
        reads = [(c, s, e, strand) for c, s, e in fragments for strand in "+-"]
    else:
        reads = random_reads(n, seed)
        path = write_tagalign(tmp_path / "reads.tagAlign", reads)
    return path, reads


def gzip_file(path, suffix=".gz"):
    with open(path, "rb") as src, gzip.open(path + suffix, "wb") as dst:
        shutil.copyfileobj(src, dst)
    return path + suffix


def truncated_copy(path, fraction=0.5):
    data = open(path, "rb").read()
    out = path + ".truncated.gz"
    with open(out, "wb") as f:
        f.write(data[:int(len(data) * fraction)])
    return out


def needs_tools(*tools):
    missing = [t for t in tools if shutil.which(t) is None]
    if missing:
        pytest.skip("missing tools: " + ", ".join(missing))


def read_stream(stream_fn, path):
    p = stream_fn(path)
    out = p.stdout.read()
    p.stdout.close()
    auto_shift_detect.check_returncode(p)
    return out


STREAMS = {"fragment": auto_shift_detect.fragment_to_tagalign_stream, "tagalign": auto_shift_detect.tagalign_stream}


def test_sample_reads_is_seeded_and_filters_unknown_chroms(tmp_path, genome):
    fasta, _ = genome
    reads = random_reads(400) + [("chrEBV", 10, 60, "+")]
    tagalign = write_tagalign(tmp_path / "reads.tagAlign", reads)

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        plus1, minus1 = auto_shift_detect.sample_reads(None, None, tagalign, 50, fasta, seed=7)
    assert any("chromosomes not in the reference" in str(w.message) for w in caught)
    plus2, minus2 = auto_shift_detect.sample_reads(None, None, tagalign, 50, fasta, seed=7)
    plus3, minus3 = auto_shift_detect.sample_reads(None, None, tagalign, 50, fasta, seed=8)

    assert len(plus1) + len(minus1) == 100
    assert plus1.equals(plus2) and minus1.equals(minus2)
    assert not (plus1.equals(plus3) and minus1.equals(minus3))
    sampled = set(map(tuple, plus1.values.tolist() + minus1.values.tolist()))
    assert len(sampled) == 100  # without replacement (all synthetic reads are distinct)
    known = {(c, str(s), str(e)) for c, s, e, _ in reads}
    assert sampled <= known
    assert "chrEBV" not in set(plus1["chr"]) | set(minus1["chr"])


def test_sample_reads_returns_everything_when_input_is_small(tmp_path, genome):
    fasta, _ = genome
    reads = random_reads(30)
    tagalign = write_tagalign(tmp_path / "reads.tagAlign", reads)
    plus, minus = auto_shift_detect.sample_reads(None, None, tagalign, 50, fasta, seed=1)
    assert len(plus) + len(minus) == 30
    assert len(plus) == sum(r[3] == "+" for r in reads)


def test_reservoir_sample_is_uniform():
    # every line should be picked with probability k/n
    n, k, trials = 50, 5, 4000
    lines = ["chr1\t{}\t{}\tN\t0\t+\n".format(i, i + 1).encode() for i in range(n)]

    class Stream:
        def __init__(self):
            import io
            self.stdout = io.BytesIO(b"".join(lines))

    class FakeFasta:
        def __init__(self, *a, **k):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

        def keys(self):
            return ["chr1"]

    orig = auto_shift_detect.pyfaidx.Fasta
    auto_shift_detect.pyfaidx.Fasta = FakeFasta
    try:
        counts = np.zeros(n)
        for seed in range(trials):
            sample, unknown = auto_shift_detect.sample_filtered_tagaligns(Stream(), "unused.fa", k, seed)
            assert len(sample) == k and len(set(sample)) == k and not unknown
            for line in sample:
                counts[int(line.split(b"\t")[1])] += 1
    finally:
        auto_shift_detect.pyfaidx.Fasta = orig
    expected = trials * k / n
    # binomial sd ~ sqrt(400*0.9) = 19; allow 5 sd
    assert np.all(np.abs(counts - expected) < 5 * np.sqrt(expected * (1 - k / n)))


def test_sample_reads_without_known_chroms_raises(tmp_path, genome):
    fasta, _ = genome
    tagalign = write_tagalign(tmp_path / "reads.tagAlign", [("chrUn", 10, 60, "+")])
    with pytest.raises(ValueError, match="No reads"):
        auto_shift_detect.sample_reads(None, None, tagalign, 50, fasta)


def test_check_returncode_raises_on_failure():
    auto_shift_detect.check_returncode(subprocess.Popen(["true"]))
    with pytest.raises(subprocess.CalledProcessError):
        auto_shift_detect.check_returncode(subprocess.Popen(["false"]))


def test_read_stream_reports_any_failed_process():
    # the first failure counts, wherever it is in the chain (a shell pipeline only reports its last command)
    ok = auto_shift_detect.ReadStream([subprocess.Popen(["true"]), subprocess.Popen(["true"])])
    auto_shift_detect.check_returncode(ok)
    first_failed = auto_shift_detect.ReadStream([subprocess.Popen(["false"]), subprocess.Popen(["true"])])
    with pytest.raises(subprocess.CalledProcessError) as err:
        auto_shift_detect.check_returncode(first_failed)
    assert err.value.cmd == "false | true"


@pytest.mark.parametrize("kind", ["fragment", "tagalign"])
def test_gzipped_input_streams_the_same_bytes(tmp_path, kind):
    needs_tools("awk", "gzip")
    plain, reads = write_input(tmp_path, kind, 500)
    out_plain = read_stream(STREAMS[kind], plain)
    assert read_stream(STREAMS[kind], gzip_file(plain)) == out_plain
    # a bgzip-style file (several gzip members, as from bgzip / tabix) decompresses the same way
    data = open(plain, "rb").read()
    multi = plain + ".multi.gz"
    with open(multi, "wb") as f:
        half = len(data) // 2
        f.write(gzip.compress(data[:half]) + gzip.compress(data[half:]))
    assert read_stream(STREAMS[kind], multi) == out_plain
    expected = "".join("{}\t{}\t{}\t{}\t{}\t{}\n".format(c, s, e, "1000" if kind == "fragment" else "N",
                                                           "0" if kind == "fragment" else "1000", strand)
                       for c, s, e, strand in reads)
    assert out_plain.decode() == expected


def test_fragment_stream_matches_the_chrombpnet_1x_shell_pipeline(tmp_path):
    # byte-for-byte what `zcat/cat file | awk ...` produced, including lines with too few or extra columns
    needs_tools("awk", "gzip", "sh")
    frag = str(tmp_path / "odd.tsv")
    with open(frag, "w") as f:
        f.write("chr1\t10\t60\tBC1\t3\nchr2\t5\t105\nchr1  7   77  x\nchr1\t8\n\nchr2\t1\t2\tBC\t1\textra\n")
    legacy = subprocess.run("cat " + frag + """ | awk -v OFS="\\t" '{print $1,$2,$3,1000,0,"+"; """
                            """print $1,$2,$3,1000,0,"-"}'""", shell=True, check=True, stdout=subprocess.PIPE).stdout
    assert read_stream(auto_shift_detect.fragment_to_tagalign_stream, frag) == legacy
    assert read_stream(auto_shift_detect.fragment_to_tagalign_stream, gzip_file(frag)) == legacy
    for path in (frag, frag + ".gz"):  # a failure message names the input file
        p = auto_shift_detect.fragment_to_tagalign_stream(path)
        p.stdout.read()
        p.stdout.close()
        assert p.wait() == 0 and path in p.args


@pytest.mark.parametrize("kind", ["fragment", "tagalign"])
def test_gzip_trailing_garbage_is_only_a_warning(tmp_path, kind):
    # gzip exits 2 ("trailing garbage ignored") after writing the whole file
    needs_tools("awk", "gzip")
    plain, _ = write_input(tmp_path, kind, 50)
    gz = gzip_file(plain)
    with open(gz, "ab") as f:
        f.write(b"\0" * 64 + b"not gzip")
    assert read_stream(STREAMS[kind], gz) == read_stream(STREAMS[kind], plain)


@pytest.mark.parametrize("kind", ["fragment", "tagalign"])
def test_truncated_gz_input_raises(tmp_path, genome, kind):
    # `zcat file | awk` exited 0 on a truncated fragment file, so the shift and the bigwig used part of the reads
    needs_tools("awk", "gzip")
    fasta, _ = genome
    plain, _ = write_input(tmp_path, kind, 3000)
    truncated = truncated_copy(gzip_file(plain))
    with pytest.raises(subprocess.CalledProcessError):
        read_stream(STREAMS[kind], truncated)
    inputs = {"fragment": (None, truncated, None), "tagalign": (None, None, truncated)}[kind]
    with pytest.raises(subprocess.CalledProcessError):
        auto_shift_detect.sample_reads(*inputs, 50, fasta)


@pytest.mark.needs_cli
def test_missing_bam_raises(tmp_path, genome):
    if shutil.which("bedtools") is None:
        pytest.skip("bedtools not on PATH")
    fasta, _ = genome
    with pytest.raises(subprocess.CalledProcessError):
        auto_shift_detect.sample_reads(str(tmp_path / "missing.bam"), None, None, 10, fasta)


def embed(pwm, offset, length=40):
    out = np.full((length, 4), 0.25)
    out[offset:offset + len(pwm)] = pwm / pwm.sum(-1, keepdims=True)
    return out


@pytest.mark.parametrize("shift", [(0, 0), (4, -5), (3, -6), (5, -4)])
def test_compute_shift_ATAC_on_embedded_reference(shift):
    plus, minus = auto_shift_detect.get_ref_pwms(get_default_data_path(DefaultDataFile.atac_ref_motifs))
    assert len(plus) == 3 and len(minus) == 3
    plus_pwm = embed(plus["GSE101074_naive_hESC_ATAC_plus"], 14 - shift[0])
    minus_pwm = embed(minus["GSE101074_naive_hESC_ATAC_minus"], 5 - shift[1])
    assert auto_shift_detect.compute_shift_ATAC(plus, minus, plus_pwm, minus_pwm) == shift


def test_compute_shift_ATAC_rejects_nonstandard_shift():
    plus, minus = auto_shift_detect.get_ref_pwms(get_default_data_path(DefaultDataFile.atac_ref_motifs))
    plus_pwm = embed(plus["GSE101074_naive_hESC_ATAC_plus"], 14)
    minus_pwm = embed(minus["GSE101074_naive_hESC_ATAC_minus"], 4)
    with pytest.raises(ValueError, match="non-standard"):
        auto_shift_detect.compute_shift_ATAC(plus, minus, plus_pwm, minus_pwm)


@pytest.mark.parametrize("shift", [(0, 0), (0, 1)])
def test_compute_shift_DNASE_on_embedded_reference(shift):
    plus, minus = auto_shift_detect.get_ref_pwms(get_default_data_path(DefaultDataFile.dnase_ref_motifs))
    assert len(plus) == 6 and len(minus) == 6
    plus_pwm = embed(plus["ENCSR000EMN_DNASE_plus"], 10 - shift[0])
    minus_pwm = embed(minus["ENCSR000EMN_DNASE_minus"], 10 - shift[1])
    assert auto_shift_detect.compute_shift_DNASE(plus, minus, plus_pwm, minus_pwm) == shift


def needs_bigwig_tools():
    needs_tools("bedtools", "bedGraphToBigWig", "awk", "sort", "gzip")


def expected_insertions(reads, plus_delta, minus_delta):
    # bedtools genomecov -5: + strand counts `start`, - strand counts `end - 1`
    counts = {}
    for chrom, start, end, strand in reads:
        pos = start + plus_delta if strand == "+" else end + minus_delta - 1
        counts[(chrom, pos)] = counts.get((chrom, pos), 0) + 1
    return counts


@pytest.mark.needs_cli
@pytest.mark.parametrize("no_st", [False, True])
def test_reads_to_bigwig_with_explicit_shift(tmp_path, genome, no_st):
    needs_bigwig_tools()
    fasta, sizes = genome
    reads = random_reads(300, seed=3)
    tagalign = write_tagalign(tmp_path / "reads.tagAlign", reads)
    args = argparse.Namespace(genome=fasta, input_bam_file=None, input_fragment_file=None,
                              input_tagalign_file=tagalign, chrom_sizes=sizes,
                              output_prefix=str(tmp_path / "out"), data_type="ATAC", bsort=False,
                              no_st=no_st, tmpdir=None, plus_shift=0, minus_shift=0,
                              ATAC_ref_path=None, DNASE_ref_path=None, num_samples=100, shift_seed=1234)
    reads_to_bigwig.main(args)
    assert_bigwig_counts(str(tmp_path / "out_unstranded.bw"), expected_insertions(reads, 4, -4))


def assert_bigwig_counts(path, expected):
    bw = pyBigWig.open(path)
    for chrom in ("chr1", "chr2"):
        vals = np.nan_to_num(np.array(bw.values(chrom, 0, CHROM_LEN)))
        want = np.zeros(CHROM_LEN)
        for (c, pos), n in expected.items():
            if c == chrom:
                want[pos] = n
        np.testing.assert_array_equal(vals, want)
    bw.close()


def bigwig_inputs(kind, path):
    # (input_bam_file, input_fragment_file, input_tagalign_file)
    return {"fragment": (None, path, None), "tagalign": (None, None, path)}[kind]


@pytest.mark.needs_cli
@pytest.mark.parametrize("no_st", [False, True])
@pytest.mark.parametrize("kind", ["fragment", "tagalign"])
def test_reads_to_bigwig_gz_input_matches_plain(tmp_path, genome, kind, no_st):
    needs_bigwig_tools()
    fasta, sizes = genome
    plain, reads = write_input(tmp_path, kind, 300)
    outputs = []
    for name, path in (("plain", plain), ("gz", gzip_file(plain))):
        prefix = str(tmp_path / name)
        reads_to_bigwig.generate_bigwig(*bigwig_inputs(kind, path), prefix, fasta, False, None, no_st, sizes, 4, -4)
        outputs.append(prefix + "_unstranded.bw")
    assert open(outputs[0], "rb").read() == open(outputs[1], "rb").read()
    assert_bigwig_counts(outputs[1], expected_insertions(reads, 4, -4))


@pytest.mark.needs_cli
@pytest.mark.parametrize("no_st", [False, True])
@pytest.mark.parametrize("kind", ["fragment", "tagalign"])
def test_reads_to_bigwig_truncated_gz_raises(tmp_path, genome, kind, no_st):
    needs_bigwig_tools()
    fasta, sizes = genome
    plain, _ = write_input(tmp_path, kind, 3000)
    truncated = truncated_copy(gzip_file(plain))
    with pytest.raises(subprocess.CalledProcessError, match="gzip"):
        reads_to_bigwig.generate_bigwig(*bigwig_inputs(kind, truncated), str(tmp_path / "out"), fasta, False, None,
                                        no_st, sizes, 4, -4)
    assert not (tmp_path / "out_unstranded.bw").exists()


@pytest.mark.needs_cli
@pytest.mark.parametrize("no_st", [False, True])
def test_reads_to_bigwig_fails_on_contig_missing_from_chrom_sizes(tmp_path, genome, no_st, capsys):
    # reads on chr2 but chrom sizes without it (e.g. a BAM with *_random contigs and a main-chromosome chrom sizes
    # file): bedtools genomecov only warns, so bedGraphToBigWig must fail the run with a hint and write no bigwig.
    needs_bigwig_tools()
    fasta, sizes = genome
    tagalign = write_tagalign(tmp_path / "reads.tagAlign", random_reads(200))
    only_chr1 = tmp_path / "chr1.chrom.sizes"
    only_chr1.write_text("chr1\t{}\n".format(CHROM_LEN))
    with pytest.raises(subprocess.CalledProcessError):
        reads_to_bigwig.generate_bigwig(None, None, tagalign, str(tmp_path / "out"), fasta, False, None, no_st,
                                        str(only_chr1), 4, -4)
    assert "must be listed in the chrom sizes file" in capsys.readouterr().err
    assert not (tmp_path / "out_unstranded.bw").exists()


@pytest.mark.needs_cli
def test_reads_to_bigwig_missing_chrom_sizes_raises(tmp_path, genome):
    needs_bigwig_tools()
    fasta, sizes = genome
    tagalign = write_tagalign(tmp_path / "reads.tagAlign", random_reads(20))
    with pytest.raises(subprocess.CalledProcessError):
        reads_to_bigwig.generate_bigwig(None, None, tagalign, str(tmp_path / "out"), fasta, False, None, False,
                                        str(tmp_path / "missing.chrom.sizes"), 4, -4)


@pytest.fixture
def recorded_shift_seeds(monkeypatch):
    seeds = []

    def compute_shift(*args):
        seeds.append(args[-1])
        return 4, -5

    monkeypatch.setattr(auto_shift_detect, "compute_shift", compute_shift)
    monkeypatch.setattr(reads_to_bigwig, "generate_bigwig", lambda *args: None)
    return seeds


PIPELINE_ARGV = ["pipeline", "-g", "g.fa", "-c", "c.sizes", "-ifrag", "f.tsv.gz", "-o", "out", "-d", "ATAC",
                 "-p", "p.bed", "-n", "n.bed", "-fl", "fold.json", "-b", "bias.h5"]


@pytest.mark.parametrize("argv", [PIPELINE_ARGV, ["bias"] + PIPELINE_ARGV[:-2] + ["-b", "0.5"]])
def test_pipeline_training_seed_does_not_change_the_shift_sample(recorded_shift_seeds, argv):
    # the pipelines hand their parsed args to reads_to_bigwig.main; their --seed is the training seed
    import chrombpnet.parsers as parsers
    for seed in ("1234", "99"):
        args = parsers.read_parser(argv + ["-s", seed])
        args.output_prefix, args.plus_shift, args.minus_shift = "out/auxiliary/data", None, None
        reads_to_bigwig.main(args)
    assert recorded_shift_seeds == [1234, 1234]


def test_standalone_seed_flag_sets_the_shift_sample_seed(recorded_shift_seeds, monkeypatch):
    base = ["reads_to_bigwig", "-g", "g.fa", "-ifrag", "f.tsv", "-c", "c.sizes", "-op", "out", "-d", "ATAC"]
    monkeypatch.setattr(sys, "argv", base)
    assert reads_to_bigwig.parse_args().shift_seed == 1234
    for flag in ("-s", "--seed"):
        monkeypatch.setattr(sys, "argv", base + [flag, "7"])
        args = reads_to_bigwig.parse_args()
        assert args.shift_seed == 7 and not hasattr(args, "seed")
        reads_to_bigwig.main(args)
    assert recorded_shift_seeds == [7, 7]

    base = ["auto_shift_detect", "-g", "g.fa", "-itag", "r.tagAlign", "-d", "DNASE"]
    monkeypatch.setattr(sys, "argv", base)
    assert auto_shift_detect.parse_args().shift_seed == 1234
    monkeypatch.setattr(sys, "argv", base + ["-s", "5"])
    auto_shift_detect.main()
    assert recorded_shift_seeds == [7, 7, 5]


def test_build_pwm_from_bigwig_plot(tmp_path, genome):
    from chrombpnet.helpers.preprocessing.analysis import build_pwm_from_bigwig
    import matplotlib
    matplotlib.use("Agg")

    fasta, sizes = genome
    bw_path = str(tmp_path / "cov.bw")
    bw = pyBigWig.open(bw_path, "w")
    bw.addHeader([("chr1", CHROM_LEN), ("chr2", CHROM_LEN)])
    starts = list(range(100, CHROM_LEN - 100, 37))
    bw.addEntries(["chr1"] * len(starts), starts, ends=[s + 1 for s in starts], values=[1.0 + (s % 3) for s in starts])
    bw.close()

    args = argparse.Namespace(bigwig=bw_path, genome=fasta, output_prefix=str(tmp_path / "qc"), chr="chr1",
                              chrom_sizes=sizes, pwm_width=24)
    build_pwm_from_bigwig.main(args)
    assert (tmp_path / "qc.png").stat().st_size > 0
