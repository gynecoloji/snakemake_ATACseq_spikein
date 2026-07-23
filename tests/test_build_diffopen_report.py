import importlib.util
import pathlib


def _load():
    p = pathlib.Path(__file__).resolve().parents[1] / "workflow" / "scripts" / "build_diffopen_report.py"
    spec = importlib.util.spec_from_file_location("bdr", p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_rnastable_in_modes():
    modes = _load().MODES
    assert "rnastable" in modes


def test_rnastable_ordered_before_anchor_shape():
    modes = list(_load().MODES)
    assert modes.index("rnastable") < modes.index("anchor_shape")
    assert modes.index("ctcf") < modes.index("rnastable")
