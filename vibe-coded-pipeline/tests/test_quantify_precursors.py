"""Tests for quantify_precursors.py (run with: python -m pytest tests)."""

import csv
import gzip
import hashlib
import json
import os
import shlex
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

import quantify_precursors as qp  # noqa: E402

FAKE_FLASHLFQ = os.path.join(HERE, "fake_flashlfq.py")
REAL_MZML_DIR = os.path.join(ROOT, "PXD009327", "MZML")


# ---------------------------------------------------------------------------
# synthetic mzML fixtures
# ---------------------------------------------------------------------------

def cv(acc, name, value="", unit_acc=None, unit_name=None):
    unit = f' unitAccession="{unit_acc}" unitName="{unit_name}" unitCvRef="UO"' if unit_acc else ""
    return f'<cvParam cvRef="MS" accession="{acc}" name="{name}" value="{value}"{unit}/>'


def ms1(index, scan, rt):
    return f"""
<spectrum id="controllerType=0 controllerNumber=1 scan={scan}" index="{index}" defaultArrayLength="0">
  {cv("MS:1000511", "ms level", 1)}
  <scanList count="1"><scan>{cv("MS:1000016", "scan start time", rt, "UO:0000031", "minute")}</scan></scanList>
</spectrum>"""


def ms2(index, scan, rt, mz=500.25, charge=2, activation=("MS:1000422",), title=None,
        rt_unit=("UO:0000031", "minute"), n_precursors=1, n_ions=1, charge_param=True):
    title_xml = cv("MS:1000796", "spectrum title", title) if title else ""
    names = {"MS:1000422": "beam-type collision-induced dissociation",
             "MS:1000133": "collision-induced dissociation",
             "MS:1000598": "electron transfer dissociation"}
    act = "".join(cv(a, names[a]) for a in activation)
    ion_params = cv("MS:1000744", "selected ion m/z", mz) if mz is not None else ""
    if charge_param and charge is not None:
        ion_params += cv("MS:1000041", "charge state", charge)
    ion_params += cv("MS:1000042", "peak intensity", 12345.5)
    ions = "".join(f"<selectedIon>{ion_params}</selectedIon>" for _ in range(n_ions))
    precursor = f"""<precursor spectrumRef="controllerType=0 controllerNumber=1 scan={scan - 1}">
      <isolationWindow>{cv("MS:1000827", "isolation window target m/z", mz or 0)}</isolationWindow>
      <selectedIonList count="{n_ions}">{ions}</selectedIonList>
      <activation>{cv("MS:1000045", "collision energy", 35)}{act}</activation></precursor>"""
    return f"""
<spectrum id="controllerType=0 controllerNumber=1 scan={scan}" index="{index}" defaultArrayLength="0">
  {cv("MS:1000511", "ms level", 2)}{title_xml}
  <scanList count="1"><scan>{cv("MS:1000016", "scan start time", rt, *rt_unit)}
    <userParam name="[Thermo Trailer Extra]Monoisotopic M/Z:" value="{mz}"/></scan></scanList>
  <precursorList count="{n_precursors}">{precursor * n_precursors}</precursorList>
</spectrum>"""


def mzml(run_id, spectra):
    return f"""<?xml version="1.0" encoding="utf-8"?>
<mzML xmlns="http://psi.hupo.org/ms/mzml" version="1.1.0">
 <run id="{run_id}"><spectrumList count="{len(spectra)}">{''.join(spectra)}</spectrumList></run>
</mzML>"""


def write_mzml(path, run_id, spectra, compress=False):
    data = mzml(run_id, spectra).encode()
    if compress:
        with gzip.open(path, "wb") as handle:
            handle.write(data)
    else:
        with open(path, "wb") as handle:
            handle.write(data)


@pytest.fixture
def two_runs(tmp_path):
    """Two small runs (one gzip, one plain, one with a space in its name)."""
    indir = tmp_path / "in put"
    indir.mkdir()
    write_mzml(indir / "runA.mzML.gz", "runA", [
        ms1(0, 1, 1.0),
        ms2(1, 2, 1.01, mz=500.25, charge=2),
        ms2(2, 3, 1.02, mz=600.30, charge=3),
        ms2(3, 4, 1.03, mz=700.35, charge=2, activation=("MS:1000133",)),   # CID: excluded
        ms2(4, 5, 1.04, mz=800.40, charge=None, charge_param=False),        # missing charge
        ms2(5, 6, 1.05, mz=900.45, charge=2),
    ], compress=True)
    write_mzml(indir / "run B.mzML", "run B", [
        ms1(0, 1, 60.0),
        ms2(1, 2, 61.2, mz=450.20, charge=2, rt_unit=("UO:0000010", "second")),
        ms2(2, 3, 62.4, mz=550.20, charge=2, title="custom title 3", rt_unit=("UO:0000010", "second")),
        ms2(3, 4, 63.6, mz=650.20, charge=4, rt_unit=("UO:0000010", "second")),
    ])
    return indir


# ---------------------------------------------------------------------------
# precursor mass and identifiers
# ---------------------------------------------------------------------------

def test_neutral_mass_matches_mzlib_convention():
    mz, z = 419.315765380859, 2
    mass = qp.neutral_mass(mz, z)
    assert mass == pytest.approx(2 * (mz - 1.007276466879), abs=1e-12)
    # FlashLFQ's ToMz(mass, z) = mass/z + proton must give back the selected m/z
    assert mass / z + qp.PROTON_MASS == pytest.approx(mz, abs=1e-9)
    with pytest.raises(ValueError):
        qp.neutral_mass(500.0, 0)


def test_record_id_is_file_hash_plus_scan_and_base_sequence_is_not_a_peptide():
    file_hash = hashlib.sha1(b"HUVEC_IL1b2b_sn").hexdigest()[:12]   # run name, no extension
    for source in ("HUVEC_IL1b2b_sn.mzML.gz", "HUVEC_IL1b2b_sn.mzML"):
        assert qp.synthetic_record_id(source, "2006", 5310) == f"{file_hash}_SCAN2006"
    assert qp.synthetic_record_id("HUVEC_IL1b2b_sn.mzML", "", 5310) == f"{file_hash}_INDEX5310"
    assert qp.RECORD_ID_RE.match(f"{file_hash}_SCAN2006")
    assert not qp.RECORD_ID_RE.match(f"{file_hash}_SCAN2006|{file_hash}_SCAN7")
    amino_acids = set("ACDEFGHIKLMNOPQRSTUVWY")
    assert set(qp.BASE_SEQUENCE_PLACEHOLDER) - amino_acids, \
        "Base Sequence must contain non-residue characters (forces the averagine isotope model)"


def _id_record(source, scan, index):
    return qp.PrecursorRecord(
        internal_id=f"{source}::{index}", source_file=source, run_id=source, spectrum_index=index,
        native_scan_id=f"scan={scan}", scan_number=scan, spectrum_title=f"{source}:{index}",
        spectrum_title_source="mzml", ms_level=2, precursor_mz=500.0, precursor_charge=2,
        retention_time=1.0, retention_time_original=1.0, retention_time_original_unit="minute",
        activation_method="HCD", hcd_evidence="", precursor_intensity=None, isolation_target_mz=None,
        thermo_monoisotopic_mz=None, n_precursors=1, n_selected_ions=1, precursor_spectrum_ref="",
    )


def test_record_ids_per_file_scan_or_index_fallback():
    records = [_id_record("a.mzML", "7", 0), _id_record("a.mzML", "9", 1),
               _id_record("b.mzML", "7", 0), _id_record("b.mzML", "7", 1),   # repeated scan: use index
               _id_record("c.mzML", "", 0)]                                  # no scan number: use index
    qp.assign_synthetic_ids(records)
    h = {name: hashlib.sha1(name.encode()).hexdigest()[:12] for name in "abc"}
    assert [r.synthetic_record_id for r in records] == [
        f"{h['a']}_SCAN7", f"{h['a']}_SCAN9", f"{h['b']}_INDEX0", f"{h['b']}_INDEX1", f"{h['c']}_INDEX0"]


def test_number_format_has_no_exponent():
    assert qp.format_number(1e-7) == "0.000000"
    assert "e" not in qp.format_number(123456789.123456789).lower()


# ---------------------------------------------------------------------------
# extraction
# ---------------------------------------------------------------------------

def test_extraction_selects_hcd_and_preserves_metadata(two_runs):
    a = qp.extract_file(str(two_runs / "runA.mzML.gz"), "{run_id}: {native_id}")
    assert a.spectrum_counts["hcd_ms2"] == 4
    assert a.spectrum_counts["ms2_not_hcd"] == 1
    assert a.spectrum_counts["ms_level_1"] == 1
    first = a.records[0]
    assert first.spectrum_title == "runA: controllerType=0 controllerNumber=1 scan=2"
    assert first.spectrum_title_source == "reconstructed"
    assert first.internal_id == "runA.mzML.gz::controllerType=0 controllerNumber=1 scan=2"
    assert first.scan_number == "2"
    assert first.precursor_mz == 500.25 and first.precursor_charge == 2
    assert first.retention_time == pytest.approx(1.01)
    assert first.precursor_intensity == 12345.5
    assert "beam-type" in first.activation_method
    statuses = {r.scan_number: r.precursor_status for r in a.records}
    assert statuses == {"2": "ok", "3": "ok", "5": "missing_precursor_charge", "6": "ok"}


def test_extraction_converts_seconds_and_keeps_mzml_title(two_runs):
    b = qp.extract_file(str(two_runs / "run B.mzML"), "{run_id}: {native_id}")
    rec = {r.scan_number: r for r in b.records}
    assert rec["2"].retention_time == pytest.approx(61.2 / 60.0)
    assert rec["2"].retention_time_original == 61.2
    assert rec["2"].retention_time_original_unit == "second"
    assert rec["3"].spectrum_title == "custom title 3"
    assert rec["3"].spectrum_title_source == "mzml_spectrum_title"


def test_ambiguous_precursors_are_flagged_not_guessed(tmp_path):
    path = tmp_path / "amb.mzML"
    write_mzml(path, "amb", [
        ms2(0, 1, 1.0, n_precursors=2),
        ms2(1, 2, 1.1, n_ions=2),
        ms2(2, 3, 1.2, mz=None),
        ms2(3, 4, 1.3, activation=("MS:1000422", "MS:1000598")),  # EThcD-like: excluded
    ])
    result = qp.extract_file(str(path), "{run_id}: {native_id}")
    statuses = [r.precursor_status for r in result.records]
    assert statuses == ["ambiguous_multiple_precursors", "ambiguous_multiple_selected_ions", "missing_precursor_mz"]
    assert result.spectrum_counts["ms2_not_hcd"] == 1
    assert all(r.neutral_mass is None for r in result.records)


def test_duplicate_native_ids_abort(tmp_path):
    path = tmp_path / "dup.mzML"
    write_mzml(path, "dup", [ms2(0, 7, 1.0), ms2(1, 7, 1.1)])
    with pytest.raises(qp.PipelineError, match="Duplicate native spectrum id"):
        qp.extract_file(str(path), "{run_id}: {native_id}")


def test_discovery_rejects_same_run_twice(tmp_path):
    write_mzml(tmp_path / "x.mzML", "x", [])
    write_mzml(tmp_path / "x.mzML.gz", "x", [], compress=True)
    with pytest.raises(qp.PipelineError, match="share a run name"):
        qp.discover_mzml_files(str(tmp_path))


# ---------------------------------------------------------------------------
# FlashLFQ input and command line
# ---------------------------------------------------------------------------

def _input_file(tmp_path, rows):
    path = tmp_path / "flashlfq_input.tsv"
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(qp.FLASHLFQ_INPUT_HEADER)
        writer.writerows(rows)
    return str(path)


def test_flashlfq_input_validation(tmp_path):
    good = ["r.mzML", "1.000000", "2", qp.BASE_SEQUENCE_PLACEHOLDER, "4dc7c9ec434e_SCAN12", "998.485447",
            qp.PROTEIN_PLACEHOLDER]
    assert qp.validate_flashlfq_input(_input_file(tmp_path, [good]), {"r.mzML"}) == {"r.mzML": 1}
    for bad in (
        [*good[:5], "9.98e2", good[6]],                      # exponent: FlashLFQ would drop it
        ["other.mzML", *good[1:]],                          # unknown spectra file
        good[:4] + ["4dc7c9ec434e_SCAN1|4dc7c9ec434e_SCAN2"] + good[5:],    # separator inside id
        good[:4] + ["F001_SCAN12"] + good[5:],               # not <12 hex>_SCAN<n>
        good[:3] + ["PEPTIDE"] + good[4:],                  # valid residues: peptide isotope model
    ):
        with pytest.raises(qp.PipelineError):
            qp.validate_flashlfq_input(_input_file(tmp_path, [bad]), {"r.mzML"})
    with pytest.raises(qp.PipelineError, match="duplicate"):
        qp.validate_flashlfq_input(_input_file(tmp_path, [good, good]), {"r.mzML"})


def test_cli_args_never_pass_false_to_switches():
    params = qp.FlashLfqParams(ppm=10, isotope_ppm=5, num_isotopes_required=2,
                               id_specific_charge=False, integrate=False, threads=4)
    args = params.cli_args("in.tsv", "rep", "out")
    assert "false" not in [a.lower() for a in args], "FlashLFQ bool options are switches; 'false' enables them"
    assert "--mbr" not in args and "--nor" not in args and "--sil" not in args
    assert args[args.index("--pipfdr") + 1] == "0"
    assert "--chg" not in args and "--int" not in args
    params.id_specific_charge = params.integrate = True
    args = params.cli_args("in.tsv", "rep", "out")
    assert "--chg" in args and "--int" in args


def test_chg_and_int_are_on_by_default(tmp_path):
    base = ["--input", str(tmp_path), "--output", str(tmp_path / "o")]
    args = qp.parse_args(base)
    assert args.id_specific_charge and args.integrate
    args = qp.parse_args(base + ["--no-id-specific-charge", "--no-integrate"])
    assert not args.id_specific_charge and not args.integrate


def test_settings_verification():
    params = qp.FlashLfqParams(10, 5, 2, False, False, 4)
    good = dict(params.expected_settings())
    assert qp.settings_mismatches(good, params.expected_settings()) == {}
    bad = dict(good, Normalize=True, PpmTolerance=20.0)
    assert set(qp.settings_mismatches(bad, params.expected_settings())) == {"Normalize", "PpmTolerance"}


def test_docker_command_mounts_and_quoting():
    params = qp.FlashLfqParams(10, 5, 2, False, False, 16)
    cmd = qp.build_docker_command(
        "img@sha256:abc", "/tmp/stage dir/run001", "run,1.mzML", "/data/in put/run,1.mzML",
        "/data/out/flashlfq/run,1", params, 16, "8g", [],
    )
    text = shlex.join(cmd)
    assert cmd[:3] == ["docker", "run", "--rm"]
    assert "--network" in cmd and cmd[cmd.index("--network") + 1] == "none"
    mounts = [cmd[i + 1] for i, a in enumerate(cmd) if a == "--mount"]
    assert mounts[0] == "type=bind,source=/tmp/stage dir/run001,target=/data/mzml"
    assert mounts[1] == 'type=bind,"source=/data/in put/run,1.mzML","target=/data/mzml/run,1.mzML",readonly'
    assert mounts[2].endswith("target=/data/out")
    assert cmd[cmd.index("--cpus") + 1] == "16" and cmd[cmd.index("--memory") + 1] == "8g"
    assert cmd[cmd.index("img@sha256:abc") + 1:][:2] == ["--idt", f"/data/out/{qp.RUN_INPUT_NAME}"]
    assert "--privileged" not in text


# ---------------------------------------------------------------------------
# FlashLFQ output mapping
# ---------------------------------------------------------------------------

def _record(source, scan, rid, rt=10.0, charge=2, mass=1000.0):
    return qp.PrecursorRecord(
        internal_id=f"{source}::scan={scan}", source_file=source, run_id=source, spectrum_index=scan,
        native_scan_id=f"scan={scan}", scan_number=str(scan), spectrum_title=f"{source}: scan={scan}",
        spectrum_title_source="reconstructed", ms_level=2, precursor_mz=mass / charge + qp.PROTON_MASS,
        precursor_charge=charge, retention_time=rt, retention_time_original=rt,
        retention_time_original_unit="minute", activation_method="HCD", hcd_evidence="activation_cvParam",
        precursor_intensity=None, isolation_target_mz=None, thermo_monoisotopic_mz=None, n_precursors=1,
        n_selected_ions=1, precursor_spectrum_ref="", neutral_mass=mass, synthetic_record_id=rid,
    )


def _feature(file_name, ids, apex=True, intensity=100.0, psms=None):
    return qp.Feature(
        feature_id="", file_name=file_name, record_ids=ids, intensity=intensity if apex else 0.0,
        rt_start=9.9 if apex else None, rt_apex=10.0 if apex else None, rt_end=10.1 if apex else None,
        fwhm=None, fwhm_status="", peak_mz=501.0 if apex else None, peak_charge=2 if apex else None,
        n_charge_states=1, detection_type="MSMS", psms_mapped=len(ids) if psms is None else psms,
        mass_error_ppm=None, first_record_mass=None, first_record_charge=None, first_record_rt=None,
    )


def test_merged_peak_maps_back_to_every_spectrum():
    records = [_record("a.mzML", 1, "F001_SCAN1"), _record("a.mzML", 2, "F001_SCAN2"), _record("a.mzML", 3, "F001_SCAN3")]
    features = [_feature("a", ["F001_SCAN1", "F001_SCAN2"]), _feature("a", ["F001_SCAN3"], apex=False)]
    qp.assign_feature_ids(features)
    mapping, problems = qp.map_features(records, features, {"a.mzML": "a.mzML"})
    assert problems == []
    assert mapping["F001_SCAN1"] is mapping["F001_SCAN2"]
    assert mapping["F001_SCAN1"].feature_id == "a#F000001"
    assert qp.quant_status(mapping["F001_SCAN1"]) == qp.STATUS_QUANTIFIED
    assert qp.quant_status(mapping["F001_SCAN3"]) == qp.STATUS_NO_FEATURE
    assert mapping["F001_SCAN3"].feature_id == ""
    assert qp.quant_status(None) == qp.STATUS_MISSING_OUTPUT
    assert qp.quant_status(_feature("a", ["x"], intensity=0.0)) == qp.STATUS_ZERO


def test_mapping_detects_integrity_problems():
    records = [_record("a.mzML", 1, "F001_SCAN1"), _record("b.mzML", 1, "F002_SCAN1")]
    staged = {"a.mzML": "a.mzML", "b.mzML": "b.mzML"}
    _, problems = qp.map_features(records, [_feature("a", ["F002_SCAN1"])], staged)
    assert any(p.startswith("cross-file") for p in problems)
    _, problems = qp.map_features(records, [_feature("a", ["F001_SCAN1"]), _feature("a", ["F001_SCAN1"])], staged)
    assert any("more than one" in p for p in problems)
    _, problems = qp.map_features(records, [_feature("a", ["F001_SCAN1"], psms=3)], staged)
    assert any("PSMs Mapped" in p for p in problems)
    _, problems = qp.map_features(records, [_feature("a", ["F009_SCAN9"])], staged)
    assert any("unknown id" in p for p in problems)


# ---------------------------------------------------------------------------
# orchestration with a fake FlashLFQ (no Docker needed)
# ---------------------------------------------------------------------------

def _run(two_runs, tmp_path, *extra):
    out = tmp_path / "out"
    local = shlex.join([sys.executable, FAKE_FLASHLFQ])
    rc = qp.main(["--input", str(two_runs), "--output", str(out), "--threads", "1",
                  "--flashlfq-local-cmd", local, "--tmp-dir", str(tmp_path / "tmp"), *extra])
    return rc, out


def _read(path):
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_end_to_end_with_fake_flashlfq(two_runs, tmp_path):
    rc, out = _run(two_runs, tmp_path)
    summary = json.load(open(out / "quantification_summary.json"))
    assert rc == 0, {k: v for k, v in summary["checks"].items() if not v["passed"]}
    rows = _read(out / "ms1_precursor_quantification.tsv")
    inventory = _read(out / "precursor_inventory.tsv")
    assert len(rows) == len(inventory) == 7          # every HCD spectrum, CID excluded
    assert [r["spectrum_title"] for r in rows] == [r["spectrum_title"] for r in inventory]
    by_title = {r["spectrum_title"]: r for r in rows}
    a2 = by_title["runA: controllerType=0 controllerNumber=1 scan=2"]
    a3 = by_title["runA: controllerType=0 controllerNumber=1 scan=3"]
    assert a2["feature_id"] == a3["feature_id"] != ""           # merged peak shared by two spectra
    assert a2["feature_n_spectra"] == "2" and a2["ms1_intensity"] == a3["ms1_intensity"]
    assert by_title["runA: controllerType=0 controllerNumber=1 scan=5"]["quantification_status"] == \
        "not_quantified_missing_precursor_charge"
    assert by_title["runA: controllerType=0 controllerNumber=1 scan=6"]["quantification_status"] == \
        qp.STATUS_NO_FEATURE
    assert by_title["runA: controllerType=0 controllerNumber=1 scan=6"]["ms1_intensity"] == ""
    assert by_title["custom title 3"]["source_file"] == "run B.mzML"
    for r in rows:                                   # no cross-file features
        if r["feature_id"]:
            assert r["feature_id"].startswith(qp.strip_mzml_suffix(r["source_file"]) + "#")
    assert summary["counts"]["quantification_status"] == {
        "quantified": 4, "no_ms1_feature_found": 2, "not_quantified_missing_precursor_charge": 1}
    assert not os.listdir(tmp_path / "tmp"), "staging must be cleaned up"


def test_flashlfq_failure_is_reported_per_spectrum(two_runs, tmp_path, monkeypatch):
    monkeypatch.setenv("FAKE_FLASHLFQ_FAIL", "1")
    rc, out = _run(two_runs, tmp_path)
    assert rc == 2
    rows = _read(out / "ms1_precursor_quantification.tsv")
    assert len(rows) == 7
    statuses = {r["quantification_status"] for r in rows}
    assert statuses == {qp.STATUS_RUN_FAILED, "not_quantified_missing_precursor_charge"}
    summary = json.load(open(out / "quantification_summary.json"))
    assert not summary["checks"]["flashlfq_runs_succeeded"]["passed"]


def test_flashlfq_timeout_marks_run_failed(two_runs, tmp_path, monkeypatch):
    import time
    monkeypatch.setenv("FAKE_FLASHLFQ_SLEEP", "60")
    start = time.monotonic()
    rc, out = _run(two_runs, tmp_path, "--timeout", "2")
    assert time.monotonic() - start < 40
    assert rc == 2
    statuses = {r["quantification_status"] for r in _read(out / "ms1_precursor_quantification.tsv")}
    assert qp.STATUS_RUN_FAILED in statuses and qp.STATUS_QUANTIFIED not in statuses


def test_project_folder_derives_paths(tmp_path):
    project = tmp_path / "PXD000001"
    (project / "MZML").mkdir(parents=True)
    (project / "MGF").mkdir()
    args = qp.parse_args(["--project", str(project)])
    assert args.input == str(project / "MZML")
    assert args.output == str(project / "results")
    assert args.reference_mgf_dir is None                     # MGF folder without .mgf files is ignored
    (project / "MGF" / "a.mgf.gz").write_bytes(b"")
    args = qp.parse_args(["--project", str(project)])
    assert args.reference_mgf_dir == str(project / "MGF")
    args = qp.parse_args(["--project", str(project), "--no-reference-mgf", "--output", str(tmp_path / "o")])
    assert args.reference_mgf_dir is None and args.output == str(tmp_path / "o")
    with pytest.raises(SystemExit):
        qp.parse_args(["--project", str(tmp_path / "missing")])
    with pytest.raises(SystemExit):
        qp.parse_args(["--input", str(project / "MZML")])      # no --output and no --project


def test_dry_run_writes_nothing(two_runs, tmp_path, capsys):
    out = tmp_path / "out"
    assert qp.main(["--input", str(two_runs), "--output", str(out), "--dry-run"]) == 0
    assert not out.exists()
    text = capsys.readouterr().out
    assert "docker run --rm" in text and "--pipfdr 0" in text


# ---------------------------------------------------------------------------
# real data (skipped when PXD009327 is not present)
# ---------------------------------------------------------------------------

@pytest.mark.skipif(not os.path.isdir(REAL_MZML_DIR), reason="PXD009327/MZML not available")
def test_pxd009327_hcd_counts():
    expected = {
        "VEL180306-01_FWeiland180306H01_FW_180307072156.mzML.gz": 10238,
        "VEL180306-03_Fweiland180306H02_FW.mzML.gz": 9554,
        "VEL180306-05_Fweiland180306H03_FW.mzML.gz": 9616,
    }
    for name, count in expected.items():
        result = qp.extract_file(os.path.join(REAL_MZML_DIR, name), "{run_id}: {native_id}")
        assert len(result.records) == count
        assert all(r.precursor_status == "ok" for r in result.records)
        assert result.records[0].spectrum_title.startswith(result.run_id + ": controllerType=0")
