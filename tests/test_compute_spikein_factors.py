import importlib.util, pathlib
import pytest

_p = pathlib.Path(__file__).resolve().parents[1] / "workflow" / "scripts" / "compute_spikein_factors.py"
_spec = importlib.util.spec_from_file_location("compute_spikein_factors", _p)
csf = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(csf)


def test_compute_factors_min_gets_one():
    f = csf.compute_factors({"a": 100, "b": 200, "c": 400})
    assert f["a"] == 1.0
    assert f["b"] == 0.5
    assert f["c"] == 0.25


def test_compute_factors_zero_raises():
    with pytest.raises(ValueError):
        csf.compute_factors({"a": 0, "b": 10})


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
