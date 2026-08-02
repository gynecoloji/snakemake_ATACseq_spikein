import importlib.util, pathlib
import pytest

_p = pathlib.Path(__file__).resolve().parents[1] / "workflow" / "scripts" / "compute_spikein_factors.py"
_spec = importlib.util.spec_from_file_location("compute_spikein_factors", _p)
csf = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(csf)


# --- factor arithmetic (floor disabled so the ratios are the only thing tested) ---

def test_compute_factors_min_gets_one():
    f = csf.compute_factors({"a": 100, "b": 200, "c": 400}, min_reads=0)
    assert f["a"] == 1.0
    assert f["b"] == 0.5
    assert f["c"] == 0.25


def test_compute_factors_zero_raises():
    with pytest.raises(ValueError):
        csf.compute_factors({"a": 0, "b": 10}, min_reads=0)


def test_compute_factors_empty_raises():
    with pytest.raises(ValueError):
        csf.compute_factors({}, min_reads=0)


# --- the depth floor ------------------------------------------------------------

def test_shallow_spikein_raises():
    """A spike-in too thin to normalize on must be an error, not a warning."""
    with pytest.raises(ValueError, match="too shallow"):
        csf.compute_factors({"a": 1_000, "b": 500_000}, min_reads=100_000)


def test_shallow_spikein_names_the_offenders():
    with pytest.raises(ValueError) as e:
        csf.compute_factors({"good": 500_000, "bad": 20}, min_reads=100_000)
    msg = str(e.value)
    assert "bad=20" in msg
    assert "good" not in msg.split("(minimum")[0]


def test_regression_gse174272_counts_are_rejected():
    """Regression: the exact counts that silently produced 65x false positives.

    GSE174272's deposited FASTQs contain no Drosophila reads. Before the floor,
    these four counts yielded a complete set of normalization factors spanning
    1.8x -- pure Poisson noise -- and 7,488 peaks at padj<0.05 where the
    spike-in-free baseline found 115.
    """
    counts = {"A485-2h": 20, "A485-6h": 20, "DMSO-1": 36, "DMSO-2": 32}
    with pytest.raises(ValueError, match="too shallow"):
        csf.compute_factors(counts)


def test_adequate_spikein_passes():
    """Real counts from GSE148175, which has a genuine Drosophila spike-in."""
    counts = {"DMSO-1": 5_933_994, "DMSO-2": 20_984_886,
              "dTAG-1": 7_852_084, "dTAG-2": 17_534_470}
    f = csf.compute_factors(counts)
    assert f["DMSO-1"] == pytest.approx(1.0)
    assert f["DMSO-2"] == pytest.approx(5_933_994 / 20_984_886)


def test_floor_can_be_disabled():
    f = csf.compute_factors({"a": 10, "b": 20}, min_reads=0)
    assert f["a"] == 1.0


# --- IO -------------------------------------------------------------------------

def test_read_count_and_sample_from_path(tmp_path):
    c = tmp_path / "GSF-Control_1.spikein_count.txt"
    c.write_text("12345\n")
    assert csf.read_count(c) == 12345
    assert csf.sample_from_path(c) == "GSF-Control_1"


def test_write_table_roundtrip(tmp_path):
    out = tmp_path / "nf.tsv"
    csf.write_table({"a": 100, "b": 200}, {"a": 1.0, "b": 0.5}, out)
    text = out.read_text().strip().splitlines()
    assert text[0] == "sample\tspikein_reads\tnorm_factor"
    assert text[1] == "a\t100\t1.000000"
    assert text[2] == "b\t200\t0.500000"


def test_main_end_to_end(tmp_path):
    for name, n in [("s1", 500_000), ("s2", 1_000_000)]:
        (tmp_path / f"{name}.spikein_count.txt").write_text(f"{n}\n")
    out = tmp_path / "nf.tsv"
    csf.main(sorted(str(p) for p in tmp_path.glob("*.spikein_count.txt")), str(out))
    lines = out.read_text().strip().splitlines()
    assert lines[1].startswith("s1\t500000\t1.000000")
