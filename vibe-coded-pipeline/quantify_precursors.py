#!/usr/bin/env python3
"""Identification-free MS1 precursor quantification of HCD MS2 spectra with FlashLFQ.

Workflow
--------
1. Stream every mzML(.gz) file in the input directory and extract the precursor
   metadata of every HCD MS2 spectrum (precursor_inventory.tsv).
2. Turn every quantifiable precursor event into one *synthetic* FlashLFQ
   "Generic" identification record (flashlfq_input.tsv). The record carries the
   neutral precursor mass, charge and MS2 retention time. "Full Sequence"
   holds a unique record id, <SHA-1 of the run name, 12 hex>_SCAN<scan number>;
   "Base Sequence" holds the constant NO_SEQUENCE, which is not a peptide, so
   FlashLFQ models every isotope envelope with averagine.
3. Run the unmodified FlashLFQ command-line program inside Docker (MBR off,
   normalisation off, no q-value/decoy columns, so no identification-based
   filtering).
4. Map every row of FlashLFQ's QuantifiedPeaks.tsv back to the contributing
   records via the placeholder identifiers that FlashLFQ writes in the
   "Full Sequence" column ("|"-separated when FlashLFQ merged several MS2
   events into one chromatographic peak).
5. Write one output row per eligible HCD spectrum
   (ms1_precursor_quantification.tsv), plus feature_mapping.tsv,
   quantification_summary.json and logs.

See README.md for the scientific rationale and the FlashLFQ source analysis.
Filesystem paths are handled with os.path on purpose (project convention).
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import logging
import os
import platform
import re
import shlex
import shutil
import statistics
import subprocess
import sys
import tempfile
import threading
import time
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, field, fields
from typing import Iterable, Iterator

PIPELINE_VERSION = "1.0.0"

# FlashLFQ 2.1.5 (git 754dbd9859ae2ac6d133532ee507ef071e5418d1, mzLib 1.0.592).
# The digest pins the exact amd64 image published on Docker Hub for tag 2.1.5.
DEFAULT_IMAGE = (
    "smithchemwisc/flashlfq:2.1.5"
    "@sha256:5340791eba72ce6c2f5d82e8d6d43a99b6498eb4de9244c213566101448adb7a"
)

# Identical to mzLib Chemistry.Constants.ProtonMass, which FlashLFQ uses for
# every m/z <-> mass conversion.
PROTON_MASS = 1.007276466879

MZML_NS = "{http://psi.hupo.org/ms/mzml}"

# PSI-MS accessions
MS_LEVEL = "MS:1000511"
SCAN_START_TIME = "MS:1000016"
FILTER_STRING = "MS:1000512"
SPECTRUM_TITLE = "MS:1000796"
SELECTED_ION_MZ = "MS:1000744"
CHARGE_STATE = "MS:1000041"
POSSIBLE_CHARGE_STATE = "MS:1000633"
PEAK_INTENSITY = "MS:1000042"
ISOLATION_TARGET_MZ = "MS:1000827"

HCD_TERMS = {
    "MS:1000422": "beam-type collision-induced dissociation",
    "MS:1002481": "higher energy beam-type collision-induced dissociation",
}
# Dissociation methods that make a spectrum "not HCD-only" when present.
OTHER_DISSOCIATION_TERMS = {
    "MS:1000133": "collision-induced dissociation",
    "MS:1000250": "electron capture dissociation",
    "MS:1000262": "infrared multiphoton dissociation",
    "MS:1000435": "photodissociation",
    "MS:1000598": "electron transfer dissociation",
    "MS:1000599": "pulsed q dissociation",
    "MS:1002472": "trap-type collision-induced dissociation",
    "MS:1002631": "Electron-Transfer/Higher-Energy Collision Dissociation (EThcD)",
    "MS:1002678": "supplemental beam-type collision-induced dissociation",
    "MS:1002679": "supplemental collision-induced dissociation",
}
UNIT_MINUTE = {"UO:0000031", "minute"}
UNIT_SECOND = {"UO:0000010", "second"}

# Precursor statuses (inventory) and quantification statuses (final table).
PRECURSOR_OK = "ok"
STATUS_QUANTIFIED = "quantified"
STATUS_NO_FEATURE = "no_ms1_feature_found"
STATUS_ZERO = "feature_with_zero_intensity"
STATUS_MISSING_OUTPUT = "failed_missing_from_flashlfq_output"
STATUS_RUN_FAILED = "failed_flashlfq_run"

FLASHLFQ_INPUT_HEADER = [
    "File Name",
    "Scan Retention Time",
    "Precursor Charge",
    "Base Sequence",
    "Full Sequence",
    "Peptide Monoisotopic Mass",
    "Protein Accession",
]
# Constant placeholder; FlashLFQ needs the column but no protein inference is
# performed or used by this pipeline (QuantifiedProteins.tsv is ignored).
PROTEIN_PLACEHOLDER = "NO_PROTEIN_INFERENCE"
# FlashLFQ splits the two sequence columns by role (mzLib FlashLfqEngine):
# - "Base Sequence" only selects the isotope model. If every character is an
#   amino-acid residue, the envelope is computed from that peptide's formula.
#   The '_' makes this constant invalid, so the envelope is averagine scaled to
#   the supplied mass, the neutral model for a precursor of unknown composition.
# - "Full Sequence" is the identity key (isotope cache, grouping, output). It
#   holds the unique record id, which FlashLFQ echoes in QuantifiedPeaks.tsv.
BASE_SEQUENCE_PLACEHOLDER = "NO_SEQUENCE"
RECORD_ID_HASH_LENGTH = 12
RECORD_ID_RE = re.compile(r"^[0-9a-f]{12}_(SCAN|INDEX)\d+$")
FLASHLFQ_OUTPUT_FILES = (
    "QuantifiedPeaks.tsv",
    "QuantifiedPeptides.tsv",
    "QuantifiedProteins.tsv",
    "BayesianFoldChangeAnalysis.tsv",
    "FlashLfqSettings.toml",
    "ErrorReport.txt",
)
# The FlashLFQ CLI returns exit code 0 even when it aborts, so its console
# output must be checked for these messages.
FLASHLFQ_FAILURE_MARKERS = (
    "Error:",
    "Problem reading PSMs",
    "FlashLFQ has crashed",
    "Could not write FlashLFQ output",
    "No peptide IDs for the specified spectra files were found",
    "Could not find experimental design file",
    "Multiple spectra files with the same name",
    "contained no MS1 peaks",
)

CONTAINER_MZML_DIR = "/data/mzml"
CONTAINER_OUT_DIR = "/data/out"

log = logging.getLogger("ms1quant")


class PipelineError(RuntimeError):
    """Fatal, user-facing pipeline error."""


# --------------------------------------------------------------------------
# Data structures
# --------------------------------------------------------------------------


@dataclass
class PrecursorRecord:
    """One HCD MS2 spectrum and its precursor metadata (one inventory row)."""

    internal_id: str
    source_file: str
    run_id: str
    spectrum_index: int
    native_scan_id: str
    scan_number: str
    spectrum_title: str
    spectrum_title_source: str
    ms_level: int
    precursor_mz: float | None
    precursor_charge: int | None
    retention_time: float | None
    retention_time_original: float | None
    retention_time_original_unit: str
    activation_method: str
    hcd_evidence: str
    precursor_intensity: float | None
    isolation_target_mz: float | None
    thermo_monoisotopic_mz: float | None
    n_precursors: int
    n_selected_ions: int
    precursor_spectrum_ref: str
    neutral_mass: float | None = None
    precursor_status: str = PRECURSOR_OK
    synthetic_record_id: str = ""


@dataclass
class FileExtraction:
    """Result of streaming one mzML file."""

    source_file: str
    run_id: str
    title_source_counts: dict[str, int]
    spectrum_counts: dict[str, int]
    records: list[PrecursorRecord] = field(default_factory=list)


@dataclass
class Feature:
    """One row of FlashLFQ's QuantifiedPeaks.tsv."""

    feature_id: str
    file_name: str
    record_ids: list[str]
    intensity: float | None
    rt_start: float | None
    rt_apex: float | None
    rt_end: float | None
    fwhm: float | None
    fwhm_status: str
    peak_mz: float | None
    peak_charge: int | None
    n_charge_states: int | None
    detection_type: str
    psms_mapped: int | None
    mass_error_ppm: float | None
    first_record_mass: float | None
    first_record_charge: int | None
    first_record_rt: float | None


# --------------------------------------------------------------------------
# mzML extraction
# --------------------------------------------------------------------------


def open_maybe_gzip(path: str):
    """Open a file for binary reading, transparently decompressing gzip."""
    with open(path, "rb") as handle:
        magic = handle.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rb")
    return open(path, "rb")


def strip_mzml_suffix(file_name: str) -> str:
    """'run.mzML.gz' -> 'run', 'run.mzML' -> 'run' (case-insensitive)."""
    lower = file_name.lower()
    for suffix in (".mzml.gz", ".mzml"):
        if lower.endswith(suffix):
            return file_name[: -len(suffix)]
    return file_name


def discover_mzml_files(input_dir: str) -> list[str]:
    """Return sorted absolute paths of *.mzML and *.mzML.gz files."""
    if not os.path.isdir(input_dir):
        raise PipelineError(f"Input directory does not exist: {input_dir}")
    paths = []
    for name in sorted(os.listdir(input_dir)):
        full = os.path.join(input_dir, name)
        if os.path.isfile(full) and name.lower().endswith((".mzml", ".mzml.gz")):
            paths.append(os.path.abspath(full))
    if not paths:
        raise PipelineError(f"No .mzML or .mzML.gz files found in {input_dir}")
    stems = Counter(strip_mzml_suffix(os.path.basename(p)).lower() for p in paths)
    duplicates = sorted(stem for stem, n in stems.items() if n > 1)
    if duplicates:
        raise PipelineError(
            "Several input files share a run name (e.g. both X.mzML and X.mzML.gz): "
            + ", ".join(duplicates)
        )
    return paths


def _cv_params(element: ET.Element) -> dict[str, ET.Element]:
    """Direct cvParam children of an element keyed by accession."""
    return {p.get("accession"): p for p in element.findall(MZML_NS + "cvParam")}


def _float_or_none(value: str | None) -> float | None:
    if value is None or value.strip() == "":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def _scan_number(native_id: str) -> str:
    """Scan number from Thermo/Waters/Bruker-style native IDs, else ''."""
    match = re.search(r"(?:^|\s)scan=(\d+)(?:\s|$)", native_id)
    if match:
        return match.group(1)
    match = re.search(r"(?:^|\s)(?:scanId|index|spectrum)=(\d+)(?:\s|$)", native_id)
    return match.group(1) if match else ""


def _activation(precursor: ET.Element | None) -> tuple[str, bool, list[str]]:
    """Return (activation description, is_hcd_only, dissociation accessions)."""
    if precursor is None:
        return "", False, []
    activation = precursor.find(MZML_NS + "activation")
    if activation is None:
        return "", False, []
    accessions = []
    names = []
    for param in activation.findall(MZML_NS + "cvParam"):
        acc = param.get("accession", "")
        if acc in HCD_TERMS or acc in OTHER_DISSOCIATION_TERMS:
            accessions.append(acc)
            names.append(param.get("name", acc))
    is_hcd = any(a in HCD_TERMS for a in accessions) and not any(
        a in OTHER_DISSOCIATION_TERMS for a in accessions
    )
    return "+".join(names), is_hcd, accessions


def parse_spectrum(
    spectrum: ET.Element,
    source_file: str,
    run_id: str,
    title_template: str,
) -> tuple[int | None, bool, PrecursorRecord | None, str]:
    """Parse one <spectrum>.

    Returns (ms_level, is_hcd, record_or_None, reason). A record is only built
    for HCD MS2 spectra; ``reason`` explains why other spectra were skipped.
    """
    params = _cv_params(spectrum)
    ms_level_param = params.get(MS_LEVEL)
    ms_level = int(ms_level_param.get("value")) if ms_level_param is not None else None
    if ms_level != 2:
        return ms_level, False, None, "not_ms2"

    native_id = spectrum.get("id", "")
    index = int(spectrum.get("index", "-1"))
    precursors = spectrum.findall(f"{MZML_NS}precursorList/{MZML_NS}precursor")
    first_precursor = precursors[0] if precursors else None

    activation_desc, is_hcd, _ = _activation(first_precursor)
    hcd_evidence = "activation_cvParam" if is_hcd else ""

    scan = spectrum.find(f"{MZML_NS}scanList/{MZML_NS}scan")
    scan_params = _cv_params(scan) if scan is not None else {}
    filter_string = scan_params[FILTER_STRING].get("value", "") if FILTER_STRING in scan_params else ""

    if not activation_desc and "@hcd" in filter_string.lower():
        # Only used when the mzML has no activation cvParam at all.
        is_hcd = True
        activation_desc = "HCD (from filter string)"
        hcd_evidence = "filter_string"
    # Any further activation step (e.g. another precursor with ETD) disqualifies.
    for other in precursors[1:]:
        desc, other_hcd, accs = _activation(other)
        if accs and not other_hcd:
            is_hcd = False
            activation_desc = f"{activation_desc};{desc}"
    if not is_hcd:
        return ms_level, False, None, "not_hcd"

    # Retention time: scan start time of the MS2 scan, converted to minutes.
    rt_min = rt_orig = None
    rt_unit = ""
    rt_param = scan_params.get(SCAN_START_TIME)
    if rt_param is not None:
        rt_orig = _float_or_none(rt_param.get("value"))
        unit_acc = rt_param.get("unitAccession", "")
        unit_name = rt_param.get("unitName", "")
        if unit_acc in UNIT_MINUTE or unit_name in UNIT_MINUTE:
            rt_unit, rt_min = "minute", rt_orig
        elif unit_acc in UNIT_SECOND or unit_name in UNIT_SECOND:
            rt_unit = "second"
            rt_min = rt_orig / 60.0 if rt_orig is not None else None
        else:
            rt_unit = f"unknown:{unit_acc or unit_name or 'none'}"

    # Precursor / selected ion.
    selected_ions = []
    for prec in precursors:
        selected_ions.extend(prec.findall(f"{MZML_NS}selectedIonList/{MZML_NS}selectedIon"))
    mz = charge = intensity = None
    possible_charges: list[str] = []
    if selected_ions:
        ion_params = _cv_params(selected_ions[0])
        mz = _float_or_none(ion_params[SELECTED_ION_MZ].get("value")) if SELECTED_ION_MZ in ion_params else None
        if CHARGE_STATE in ion_params:
            charge_value = _float_or_none(ion_params[CHARGE_STATE].get("value"))
            charge = int(charge_value) if charge_value is not None and charge_value.is_integer() else None
        possible_charges = [
            p.get("value", "")
            for p in selected_ions[0].findall(MZML_NS + "cvParam")
            if p.get("accession") == POSSIBLE_CHARGE_STATE
        ]
        if PEAK_INTENSITY in ion_params:
            intensity = _float_or_none(ion_params[PEAK_INTENSITY].get("value"))

    isolation_target = None
    if first_precursor is not None:
        window = first_precursor.find(MZML_NS + "isolationWindow")
        if window is not None:
            wparams = _cv_params(window)
            if ISOLATION_TARGET_MZ in wparams:
                isolation_target = _float_or_none(wparams[ISOLATION_TARGET_MZ].get("value"))

    thermo_mono = None
    if scan is not None:
        for user_param in scan.findall(MZML_NS + "userParam"):
            if "Monoisotopic M/Z" in user_param.get("name", ""):
                thermo_mono = _float_or_none(user_param.get("value"))

    # Title: preserved from the mzML when present, otherwise reconstructed.
    title_param = params.get(SPECTRUM_TITLE)
    scan_number = _scan_number(native_id)
    if title_param is not None and title_param.get("value"):
        title, title_source = title_param.get("value"), "mzml_spectrum_title"
    else:
        title = title_template.format(
            run_id=run_id,
            native_id=native_id,
            file_stem=strip_mzml_suffix(source_file),
            scan_number=scan_number,
            index=index,
        )
        title_source = "reconstructed"

    # Precursor status: anything other than "ok" is reported but not quantified.
    if len(precursors) != 1:
        status = "ambiguous_multiple_precursors" if precursors else "missing_precursor"
    elif len(selected_ions) != 1:
        status = "ambiguous_multiple_selected_ions" if selected_ions else "missing_selected_ion"
    elif mz is None or mz <= 0:
        status = "missing_precursor_mz"
    elif charge is None:
        status = "ambiguous_precursor_charge" if possible_charges else "missing_precursor_charge"
    elif charge <= 0:
        status = "invalid_precursor_charge"
    elif rt_min is None or rt_min < 0:
        status = "missing_or_invalid_retention_time"
    else:
        status = PRECURSOR_OK

    record = PrecursorRecord(
        internal_id=f"{source_file}::{native_id}",
        source_file=source_file,
        run_id=run_id,
        spectrum_index=index,
        native_scan_id=native_id,
        scan_number=scan_number,
        spectrum_title=title,
        spectrum_title_source=title_source,
        ms_level=2,
        precursor_mz=mz,
        precursor_charge=charge,
        retention_time=rt_min,
        retention_time_original=rt_orig,
        retention_time_original_unit=rt_unit,
        activation_method=activation_desc,
        hcd_evidence=hcd_evidence,
        precursor_intensity=intensity,
        isolation_target_mz=isolation_target,
        thermo_monoisotopic_mz=thermo_mono,
        n_precursors=len(precursors),
        n_selected_ions=len(selected_ions),
        precursor_spectrum_ref=first_precursor.get("spectrumRef", "") if first_precursor is not None else "",
        precursor_status=status,
    )
    if status == PRECURSOR_OK:
        record.neutral_mass = neutral_mass(mz, charge)
    return ms_level, True, record, "hcd_ms2"


def neutral_mass(precursor_mz: float, charge: int) -> float:
    """Neutral mass of the selected ion: z * (m/z - proton mass).

    Equals the monoisotopic mass only if the selected m/z is the monoisotopic
    peak; FlashLFQ rejects envelopes that fit a +/-1 Da shifted mass better.
    """
    if charge <= 0:
        raise ValueError(f"charge must be positive, got {charge}")
    return charge * (precursor_mz - PROTON_MASS)


def extract_file(path: str, title_template: str) -> FileExtraction:
    """Stream one mzML(.gz) file and return all HCD MS2 precursor records.

    Uses iterparse and discards each <spectrum> after parsing so memory use is
    bounded by one spectrum, not the file size. Binary arrays are not decoded.
    """
    source_file = os.path.basename(path)
    run_id = strip_mzml_suffix(source_file)
    counts: Counter[str] = Counter()
    title_sources: Counter[str] = Counter()
    records: list[PrecursorRecord] = []
    seen_ids: set[str] = set()
    parent_stack: list[ET.Element] = []

    with open_maybe_gzip(path) as handle:
        for event, elem in ET.iterparse(handle, events=("start", "end")):
            tag = elem.tag
            if event == "start":
                if tag == MZML_NS + "run":
                    run_id = elem.get("id") or run_id
                elif tag in (MZML_NS + "spectrumList", MZML_NS + "chromatogramList"):
                    parent_stack.append(elem)
                continue
            if tag == MZML_NS + "spectrum":
                counts["spectra_total"] += 1
                native_id = elem.get("id", "")
                if native_id in seen_ids:
                    raise PipelineError(f"Duplicate native spectrum id '{native_id}' in {source_file}")
                seen_ids.add(native_id)
                ms_level, is_hcd, record, reason = parse_spectrum(elem, source_file, run_id, title_template)
                counts[f"ms_level_{ms_level}"] += 1
                if ms_level == 2 and not is_hcd:
                    counts["ms2_not_hcd"] += 1
                if record is not None:
                    counts["hcd_ms2"] += 1
                    title_sources[record.spectrum_title_source] += 1
                    records.append(record)
                if parent_stack:
                    parent_stack[-1].remove(elem)
                elem.clear()
            elif tag == MZML_NS + "chromatogram":
                if parent_stack:
                    parent_stack[-1].remove(elem)
                elem.clear()

    return FileExtraction(
        source_file=source_file,
        run_id=run_id,
        title_source_counts=dict(title_sources),
        spectrum_counts=dict(counts),
        records=records,
    )


def _extract_file_worker(args: tuple[str, str]) -> FileExtraction:
    return extract_file(*args)


def extract_all(paths: list[str], title_template: str, workers: int) -> list[FileExtraction]:
    """Extract all files, in parallel processes when workers > 1."""
    jobs = [(p, title_template) for p in paths]
    if workers <= 1 or len(paths) == 1:
        return [extract_file(*job) for job in jobs]
    with ProcessPoolExecutor(max_workers=min(workers, len(paths))) as pool:
        return list(pool.map(_extract_file_worker, jobs))


# --------------------------------------------------------------------------
# Synthetic FlashLFQ input
# --------------------------------------------------------------------------


def file_hash(source_file: str) -> str:
    """First RECORD_ID_HASH_LENGTH hex characters (lower case) of the SHA-1 of
    the run name, i.e. the file name without '.mzML' / '.mzML.gz', UTF-8 encoded.
    Reproduce with: hashlib.sha1(b"HUVEC_IL1b1a_sn").hexdigest()[:12]
    """
    run_name = strip_mzml_suffix(source_file)
    return hashlib.sha1(run_name.encode("utf-8")).hexdigest()[:RECORD_ID_HASH_LENGTH]


def synthetic_record_id(source_file: str, scan_number: str, spectrum_index: int) -> str:
    """Deterministic record id for one precursor event, written to "Full Sequence".

    '<file hash>_SCAN2006' is scan 2006 of that file, e.g. 'a1b2c3d4e5f6_SCAN2006'.
    Files without unique scan numbers use the mzML spectrum index instead:
    '<file hash>_INDEX5310'. The id contains no '|' (FlashLFQ's output
    separator) and no whitespace.
    """
    if scan_number:
        return f"{file_hash(source_file)}_SCAN{scan_number}"
    return f"{file_hash(source_file)}_INDEX{spectrum_index}"


def assign_synthetic_ids(records: list[PrecursorRecord]) -> None:
    by_file: dict[str, list[PrecursorRecord]] = defaultdict(list)
    for record in records:
        if record.precursor_status == PRECURSOR_OK:
            by_file[record.source_file].append(record)
    seen: dict[str, str] = {}
    for source, file_records in by_file.items():
        scans = [r.scan_number for r in file_records]
        use_scan = all(scans) and len(set(scans)) == len(scans)
        if not use_scan:
            log.info("event=record_ids_use_spectrum_index file=%s reason=missing or repeated scan numbers", source)
        for record in file_records:
            rid = synthetic_record_id(source, record.scan_number if use_scan else "", record.spectrum_index)
            if rid in seen:
                raise PipelineError(f"Record id collision between {seen[rid]} and {record.internal_id}")
            seen[rid] = record.internal_id
            record.synthetic_record_id = rid


def format_number(value: float, decimals: int = 6) -> str:
    """Fixed-point formatting; FlashLFQ parses with NumberStyles.Number,
    which rejects exponent notation."""
    return f"{value:.{decimals}f}"


def flashlfq_input_rows(records: Iterable[PrecursorRecord], staged_names: dict[str, str]) -> Iterator[list[str]]:
    for r in records:
        if r.precursor_status != PRECURSOR_OK:
            continue
        yield [
            staged_names[r.source_file],
            format_number(r.retention_time, 6),
            str(r.precursor_charge),
            BASE_SEQUENCE_PLACEHOLDER,
            r.synthetic_record_id,
            format_number(r.neutral_mass, 6),
            PROTEIN_PLACEHOLDER,
        ]


_NUMBER_STYLES_RE = re.compile(r"^[+-]?(\d{1,3}(,\d{3})*|\d+)(\.\d+)?$")


def validate_flashlfq_input(path: str, expected_file_names: set[str]) -> dict[str, int]:
    """Re-read flashlfq_input.tsv and apply the constraints of FlashLFQ's
    Generic PSM reader (Util/PsmReader.cs, FlashLFQ 2.1.5):

    - header recognised as "Generic" (and not MetaMorpheus/MaxQuant/...),
      without Q-Value / Target/Decoy columns, so no q-value or decoy filtering;
    - numbers parseable with .NET NumberStyles.Number (no exponent);
    - RT >= 0, integer charge > 0;
    - every Full Sequence unique (one mass per sequence; unique ids);
    - Base Sequence is the non-peptide placeholder (averagine isotope model);
    - Full Sequence is a record id without '|' or whitespace;
    - File Name maps to one staged spectra file.
    """
    file_name = os.path.basename(path)
    reserved_suffixes = ("psm.tsv", "peptide.tsv", "protein.tsv", "peaks.tsv")
    if file_name.lower().endswith(reserved_suffixes) or not file_name.lower().endswith(".tsv"):
        raise PipelineError(f"{file_name}: name would be routed to a different mzLib reader")
    seen: set[str] = set()
    per_file: Counter[str] = Counter()
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        if header != FLASHLFQ_INPUT_HEADER:
            raise PipelineError(f"Unexpected FlashLFQ input header: {header}")
        lowered = [h.lower() for h in header]
        for forbidden in ("q-value", "qvalue", "qvalue notch", "target/decoy", "decoy/contaminant/target"):
            if forbidden in lowered:
                raise PipelineError(f"FlashLFQ input must not contain column '{forbidden}'")
        for line_no, row in enumerate(reader, start=2):
            if len(row) != len(header):
                raise PipelineError(f"flashlfq_input.tsv line {line_no}: wrong column count")
            fname, rt, charge, base_seq, full_seq, mass, _protein = row
            if fname not in expected_file_names:
                raise PipelineError(f"line {line_no}: unknown spectra file '{fname}'")
            for label, value in (("RT", rt), ("charge", charge), ("mass", mass)):
                if not _NUMBER_STYLES_RE.match(value):
                    raise PipelineError(f"line {line_no}: {label} '{value}' not parseable by FlashLFQ")
            if float(rt) < 0:
                raise PipelineError(f"line {line_no}: negative retention time")
            if not float(charge).is_integer() or int(float(charge)) <= 0:
                raise PipelineError(f"line {line_no}: invalid charge {charge}")
            if base_seq != BASE_SEQUENCE_PLACEHOLDER:
                raise PipelineError(f"line {line_no}: unexpected Base Sequence {base_seq}")
            if not RECORD_ID_RE.match(full_seq):
                raise PipelineError(f"line {line_no}: unexpected record id {full_seq}")
            if full_seq in seen:
                raise PipelineError(f"line {line_no}: duplicate synthetic sequence {full_seq}")
            seen.add(full_seq)
            per_file[fname] += 1
    return dict(per_file)


# --------------------------------------------------------------------------
# FlashLFQ execution
# --------------------------------------------------------------------------


@dataclass
class FlashLfqParams:
    ppm: float
    isotope_ppm: float
    num_isotopes_required: int
    id_specific_charge: bool
    integrate: bool
    threads: int

    def cli_args(self, idt: str, rep: str, out: str) -> list[str]:
        """FlashLFQ 2.1.5 CLI arguments.

        The CLI (CommandLineParser 2.7) treats every bool option as a
        presence switch: "--nor false" ENABLES normalisation and the "false"
        token is silently ignored. Bool switches are therefore only passed
        when they must be true, and options defaulting to true (--mbr) cannot
        be switched off at all. Match-between-runs is neutralised instead:
        donors require QValue < pipfdr/5 and every synthetic record has
        QValue 0, so "--pipfdr 0" leaves no donor and no MBR peak; in addition
        each run is quantified in its own FlashLFQ invocation. The effective
        settings are verified from FlashLfqSettings.toml after every run.
        """
        args = [
            "--idt", idt,
            "--rep", rep,
            "--out", out,
            "--ppm", repr(float(self.ppm)),
            "--iso", repr(float(self.isotope_ppm)),
            "--nis", str(self.num_isotopes_required),
            "--thr", str(self.threads),
            "--pipfdr", "0",
            "--rea",
        ]
        if self.id_specific_charge:
            args.append("--chg")
        if self.integrate:
            args.append("--int")
        return args

    def expected_settings(self) -> dict[str, object]:
        """Values FlashLFQ must report in FlashLfqSettings.toml."""
        return {
            "Normalize": False,
            "PpmTolerance": float(self.ppm),
            "IsotopePpmTolerance": float(self.isotope_ppm),
            "Integrate": self.integrate,
            "NumIsotopesRequired": self.num_isotopes_required,
            "IdSpecificChargeState": self.id_specific_charge,
            "BayesianProteinQuant": False,
            "UseSharedPeptidesForProteinQuant": False,
            "MbrQValueThreshold": 0.0,
            "UsePepQValue": False,
            "RnaMode": False,
            "Silent": False,
        }


def read_flashlfq_settings(path: str) -> dict[str, object]:
    """Parse the flat key = value TOML that FlashLFQ writes."""
    settings: dict[str, object] = {}
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if "=" not in line or line.lstrip().startswith("#"):
                continue
            key, value = (part.strip() for part in line.split("=", 1))
            if value in ("true", "false"):
                settings[key] = value == "true"
            elif value.startswith('"'):
                settings[key] = value.strip('"')
            else:
                try:
                    settings[key] = int(value)
                except ValueError:
                    try:
                        settings[key] = float(value)
                    except ValueError:
                        settings[key] = value
    return settings


def settings_mismatches(actual: dict[str, object], expected: dict[str, object]) -> dict[str, object]:
    bad = {}
    for key, want in expected.items():
        got = actual.get(key)
        if isinstance(want, float) and isinstance(got, (int, float)):
            ok = abs(float(got) - want) < 1e-9
        else:
            ok = got == want
        if not ok:
            bad[key] = {"expected": want, "actual": got}
    return bad


def _mount_arg(source: str, target: str, readonly: bool) -> str:
    """--mount value, CSV-quoted so paths with spaces/commas are safe."""
    parts = ["type=bind", f"source={source}", f"target={target}"]
    if readonly:
        parts.append("readonly")
    quoted = []
    for part in parts:
        if any(ch in part for ch in ',"'):
            part = '"' + part.replace('"', '""') + '"'
        quoted.append(part)
    return ",".join(quoted)


def build_docker_command(
    image: str,
    spectra_dir: str,
    spectra_name: str,
    readonly_source: str | None,
    work_dir: str,
    params: FlashLfqParams,
    cpus: float | None,
    memory: str | None,
    extra_docker_args: list[str],
    container_name: str | None = None,
) -> list[str]:
    """docker run command for one spectra file.

    ``spectra_dir`` (a private, writable staging directory) is mounted at
    /data/mzml because FlashLFQ always writes a temporary peak index
    ("<run>.ind") next to the spectra file while MBR is enabled. An
    uncompressed input file (``readonly_source``) is bind-mounted read-only
    on top, inside that directory, so the original data cannot be modified.
    ``work_dir`` (FlashLFQ input + output) is mounted at /data/out.
    """
    cmd = [
        "docker", "run", "--rm",
        "--network", "none",
        "--security-opt", "no-new-privileges",
        "-e", "HOME=/tmp",
        "-e", "DOTNET_CLI_TELEMETRY_OPTOUT=1",
        "-e", "DOTNET_NOLOGO=1",
    ]
    if container_name:
        cmd += ["--name", container_name]
    if hasattr(os, "getuid"):
        cmd += ["--user", f"{os.getuid()}:{os.getgid()}"]
    if cpus:
        cmd += ["--cpus", str(cpus)]
    if memory:
        cmd += ["--memory", memory]
    cmd += ["--mount", _mount_arg(spectra_dir, CONTAINER_MZML_DIR, False)]
    if readonly_source:
        cmd += ["--mount", _mount_arg(readonly_source, f"{CONTAINER_MZML_DIR}/{spectra_name}", True)]
    cmd += ["--mount", _mount_arg(work_dir, CONTAINER_OUT_DIR, False)]
    cmd += extra_docker_args
    cmd.append(image)
    cmd += params.cli_args(
        idt=f"{CONTAINER_OUT_DIR}/{RUN_INPUT_NAME}",
        rep=CONTAINER_MZML_DIR,
        out=f"{CONTAINER_OUT_DIR}/{RUN_OUTPUT_DIR}",
    )
    return cmd


RUN_INPUT_NAME = "flashlfq_input_run.tsv"
RUN_OUTPUT_DIR = "output"


def is_gzip_file(path: str) -> bool:
    with open(path, "rb") as handle:
        return handle.read(2) == b"\x1f\x8b"


def stage_spectra_file(path: str, staging_dir: str, copy_plain: bool) -> tuple[str, str | None]:
    """Prepare one spectra file for FlashLFQ in its own staging directory.

    Returns (FlashLFQ-visible file name, read-only source to mount or None).
    FlashLFQ only reads uncompressed .mzML, so .mzML.gz is decompressed into
    ``staging_dir``. A plain .mzML is left in place and mounted read-only by
    Docker (returned as the read-only source), or hard-linked/copied when
    ``copy_plain`` (local development runner). Inputs are never modified.
    """
    os.makedirs(staging_dir, exist_ok=True)
    name = strip_mzml_suffix(os.path.basename(path)) + ".mzML"
    target = os.path.join(staging_dir, name)
    if is_gzip_file(path):
        log.info("event=decompress source=%s target=%s", path, target)
        with gzip.open(path, "rb") as src, open(target + ".part", "wb") as dst:
            shutil.copyfileobj(src, dst, length=16 * 1024 * 1024)
        os.replace(target + ".part", target)
        return name, None
    if copy_plain:
        try:
            os.link(path, target)
        except OSError:
            shutil.copy2(path, target)
        return name, None
    return name, path


def check_docker(image: str) -> None:
    """Fail early when Docker is unusable; pull the image if it is missing."""
    if shutil.which("docker") is None:
        raise PipelineError("'docker' executable not found on PATH")
    info = subprocess.run(["docker", "version", "--format", "{{.Server.Version}}"],
                          capture_output=True, text=True, check=False)
    if info.returncode != 0:
        raise PipelineError(f"Docker daemon not reachable: {info.stderr.strip() or info.stdout.strip()}")
    log.info("event=docker_ok server_version=%s", info.stdout.strip())
    if subprocess.run(["docker", "image", "inspect", image], capture_output=True, check=False).returncode != 0:
        log.info("event=docker_pull image=%s", image)
        pull = subprocess.run(["docker", "pull", image], capture_output=True, text=True, check=False)
        if pull.returncode != 0:
            raise PipelineError(f"docker pull {image} failed: {pull.stderr.strip()}")


def run_logged(
    cmd: list[str], log_path: str, timeout: float | None, kill_cmd: list[str] | None = None
) -> tuple[int, str]:
    """Run a command, stream stdout+stderr to a log file, return (rc, text).

    On timeout the process is killed; ``kill_cmd`` (e.g. "docker kill <name>")
    is run first, because killing the docker client alone does not stop the
    container.
    """
    log.info("event=exec command=%s", shlex.join(cmd))
    lines: list[str] = []
    start = time.monotonic()
    with open(log_path, "w", encoding="utf-8") as log_file:
        log_file.write(f"# command: {shlex.join(cmd)}\n")
        log_file.flush()
        try:
            proc = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                encoding="utf-8",
                errors="replace",
            )
        except OSError as exc:
            raise PipelineError(f"Cannot execute '{cmd[0]}': {exc}") from exc

        def pump() -> None:
            assert proc.stdout is not None
            for line in proc.stdout:
                log_file.write(line)
                log_file.flush()
                lines.append(line.rstrip("\n"))

        reader = threading.Thread(target=pump, daemon=True)
        reader.start()
        try:
            rc = proc.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            if kill_cmd:
                subprocess.run(kill_cmd, capture_output=True, check=False)
            proc.kill()
            proc.wait()
            reader.join(timeout=10)
            log_file.write(f"# killed after timeout of {timeout} s\n")
            raise PipelineError(f"FlashLFQ exceeded timeout of {timeout} s") from None
        reader.join()
        log_file.write(f"# exit code: {rc}; elapsed {time.monotonic() - start:.1f} s\n")
    return rc, "\n".join(lines)


def check_flashlfq_console(text: str, n_sent: int) -> dict[str, object]:
    """Detect FlashLFQ failures that do not change the exit code, and check
    that the reader accepted every synthetic record."""
    failures = [m for m in FLASHLFQ_FAILURE_MARKERS if m in text]
    reader_rejections = text.count("could not be read") + text.count("was not interpretable")
    found = [int(m) for m in re.findall(r"Done reading PSMs; found (\d+)", text)]
    accepted = re.findall(r"read in (\d+) identifications", text)
    version = re.search(r"FlashLFQ version: (\S+)", text)
    return {
        "failure_markers": failures,
        "reader_rejection_messages": reader_rejections,
        "records_reported_read": found[-1] if found else None,
        "records_reported_accepted": int(accepted[-1]) if accepted else None,
        "records_sent": n_sent,
        "version_line": version.group(1) if version else None,
    }


# --------------------------------------------------------------------------
# FlashLFQ output parsing and mapping
# --------------------------------------------------------------------------


def _num(value: str) -> float | None:
    value = (value or "").strip()
    if value in ("", "-", "NaN", "nan"):
        return None
    try:
        return float(value)
    except ValueError:
        return None


def _int(value: str) -> int | None:
    number = _num(value)
    return int(number) if number is not None else None


def parse_quantified_peaks(path: str) -> list[Feature]:
    """Parse QuantifiedPeaks.tsv (one row per FlashLFQ chromatographic peak).

    "Full Sequence" lists the distinct ModifiedSequences of all
    identifications merged into the peak, joined by '|'. With one unique
    synthetic sequence per MS2 event, this is the exact list of contributing
    precursor records.
    """
    features: list[Feature] = []
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"File Name", "Full Sequence", "Peak intensity", "Peak RT Apex", "PSMs Mapped"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise PipelineError(f"QuantifiedPeaks.tsv lacks columns {sorted(missing)}")
        for row in reader:
            ids = [s for s in row["Full Sequence"].split("|") if s]
            features.append(
                Feature(
                    feature_id="",
                    file_name=row["File Name"],
                    record_ids=ids,
                    intensity=_num(row["Peak intensity"]),
                    rt_start=_num(row.get("Peak RT Start", "")),
                    rt_apex=_num(row.get("Peak RT Apex", "")),
                    rt_end=_num(row.get("Peak RT End", "")),
                    fwhm=_num(row.get("Peak FWHM", "")),
                    fwhm_status=row.get("Peak FWHM Status", ""),
                    peak_mz=_num(row.get("Peak MZ", "")),
                    peak_charge=_int(row.get("Peak Charge", "")),
                    n_charge_states=_int(row.get("Num Charge States Observed", "")),
                    detection_type=row.get("Peak Detection Type", ""),
                    psms_mapped=_int(row.get("PSMs Mapped", "")),
                    mass_error_ppm=_num(row.get("Peak Apex Mass Error (ppm)", "")),
                    first_record_mass=_num(row.get("Peptide Monoisotopic Mass", "")),
                    first_record_charge=_int(row.get("Precursor Charge", "")),
                    first_record_rt=_num(row.get("MS2 Retention Time", "")),
                )
            )
    return features


def has_apex(feature: Feature) -> bool:
    return feature.rt_apex is not None


def assign_feature_ids(features: list[Feature]) -> None:
    """Deterministic ids '<run>#F000001' for peaks with an apex, numbered per
    file by apex RT, then m/z, then smallest member id. Apex-less rows are
    "no feature found" placeholders and get no feature id."""
    by_file: dict[str, list[Feature]] = defaultdict(list)
    for feature in features:
        if has_apex(feature):
            by_file[feature.file_name].append(feature)
    for file_name, items in by_file.items():
        items.sort(key=lambda f: (f.rt_apex, f.peak_mz or 0.0, min(f.record_ids)))
        for n, feature in enumerate(items, start=1):
            feature.feature_id = f"{file_name}#F{n:06d}"


def map_features(
    records: list[PrecursorRecord],
    features: list[Feature],
    staged_names: dict[str, str],
) -> tuple[dict[str, Feature], list[str]]:
    """Map synthetic record id -> FlashLFQ peak, by identifier.

    Returns (mapping, problems). Problems are integrity violations: a record
    appearing in more than one peak, unknown ids, cross-file assignments, or
    'PSMs Mapped' disagreeing with the identifier list.
    """
    record_by_id = {r.synthetic_record_id: r for r in records if r.synthetic_record_id}
    expected_run = {
        r.synthetic_record_id: strip_mzml_suffix(staged_names[r.source_file])
        for r in records
        if r.synthetic_record_id
    }
    mapping: dict[str, Feature] = {}
    problems: list[str] = []
    for feature in features:
        if feature.psms_mapped is not None and feature.psms_mapped != len(feature.record_ids):
            problems.append(
                f"peak in {feature.file_name}: PSMs Mapped={feature.psms_mapped} but "
                f"{len(feature.record_ids)} ids listed"
            )
        if len(set(feature.record_ids)) != len(feature.record_ids):
            problems.append(f"peak in {feature.file_name}: repeated ids {feature.record_ids}")
        for rid in feature.record_ids:
            if rid not in record_by_id:
                problems.append(f"unknown id {rid} in FlashLFQ output")
                continue
            if expected_run[rid] != feature.file_name:
                problems.append(f"cross-file assignment: {rid} expected {expected_run[rid]}, got {feature.file_name}")
                continue
            if rid in mapping:
                problems.append(f"record {rid} appears in more than one FlashLFQ peak")
                continue
            mapping[rid] = feature
    return mapping, problems


def quant_status(feature: Feature | None) -> str:
    if feature is None:
        return STATUS_MISSING_OUTPUT
    if not has_apex(feature):
        return STATUS_NO_FEATURE
    if not feature.intensity:
        return STATUS_ZERO
    return STATUS_QUANTIFIED


# --------------------------------------------------------------------------
# Output writing
# --------------------------------------------------------------------------


def fmt(value: object) -> str:
    """TSV cell: '' for missing, repr-precision floats, plain otherwise."""
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        return repr(value)
    return str(value)


INVENTORY_COLUMNS = [f.name for f in fields(PrecursorRecord)]

FINAL_COLUMNS = [
    # required
    "source_file", "spectrum_title", "native_scan_id", "scan_number",
    "precursor_mz", "precursor_charge", "retention_time",
    "feature_id", "feature_mz", "feature_rt_start", "feature_rt_apex", "feature_rt_end",
    "ms1_intensity", "quantification_status",
    # additional
    "internal_id", "spectrum_title_source", "precursor_status", "synthetic_record_id",
    "neutral_mass", "activation_method", "precursor_intensity", "thermo_monoisotopic_mz",
    "feature_charge", "feature_charge_matches_precursor", "feature_n_charge_states",
    "feature_n_spectra", "feature_fwhm", "feature_fwhm_status", "feature_mass_error_ppm",
    "ms2_rt_minus_feature_apex", "ms2_rt_within_feature_bounds", "flashlfq_detection_type",
]

MAPPING_COLUMNS = [
    "internal_id", "source_file", "spectrum_title", "native_scan_id", "synthetic_record_id",
    "feature_id", "feature_file", "feature_n_spectra", "feature_rt_apex", "feature_mz",
    "ms1_intensity", "quantification_status", "mapping_method",
]


def write_tsv(path: str, header: list[str], rows: Iterable[Iterable[object]]) -> int:
    tmp = path + ".part"
    n = 0
    with open(tmp, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n", quoting=csv.QUOTE_MINIMAL)
        writer.writerow(header)
        for row in rows:
            writer.writerow([fmt(v) for v in row])
            n += 1
    os.replace(tmp, path)
    return n


def final_row(record: PrecursorRecord, feature: Feature | None, status: str) -> dict[str, object]:
    row: dict[str, object] = {
        "source_file": record.source_file,
        "spectrum_title": record.spectrum_title,
        "native_scan_id": record.native_scan_id,
        "scan_number": record.scan_number,
        "precursor_mz": record.precursor_mz,
        "precursor_charge": record.precursor_charge,
        "retention_time": record.retention_time,
        "quantification_status": status,
        "internal_id": record.internal_id,
        "spectrum_title_source": record.spectrum_title_source,
        "precursor_status": record.precursor_status,
        "synthetic_record_id": record.synthetic_record_id,
        "neutral_mass": record.neutral_mass,
        "activation_method": record.activation_method,
        "precursor_intensity": record.precursor_intensity,
        "thermo_monoisotopic_mz": record.thermo_monoisotopic_mz,
    }
    if feature is not None:
        row["flashlfq_detection_type"] = feature.detection_type
    if feature is not None and has_apex(feature):
        row.update(
            feature_id=feature.feature_id,
            feature_mz=feature.peak_mz,
            feature_rt_start=feature.rt_start,
            feature_rt_apex=feature.rt_apex,
            feature_rt_end=feature.rt_end,
            ms1_intensity=feature.intensity,
            feature_charge=feature.peak_charge,
            feature_charge_matches_precursor=feature.peak_charge == record.precursor_charge,
            feature_n_charge_states=feature.n_charge_states,
            feature_n_spectra=len(feature.record_ids),
            feature_fwhm=feature.fwhm,
            feature_fwhm_status=feature.fwhm_status,
            feature_mass_error_ppm=feature.mass_error_ppm,
            ms2_rt_minus_feature_apex=record.retention_time - feature.rt_apex,
            ms2_rt_within_feature_bounds=feature.rt_start <= record.retention_time <= feature.rt_end,
        )
    return row


# --------------------------------------------------------------------------
# Validation against reference MGF titles (optional)
# --------------------------------------------------------------------------


def read_mgf_index(mgf_dir: str) -> dict[str, dict[str, str]]:
    """TITLE -> {PEPMASS, CHARGE, RTINSECONDS, file} from *.mgf(.gz)."""
    index: dict[str, dict[str, str]] = {}
    for name in sorted(os.listdir(mgf_dir)):
        if not name.lower().endswith((".mgf", ".mgf.gz")):
            continue
        path = os.path.join(mgf_dir, name)
        with open_maybe_gzip(path) as raw:
            current: dict[str, str] = {}
            for line in raw:
                text = line.decode("utf-8", "replace").rstrip("\r\n")
                if text == "BEGIN IONS":
                    current = {"file": name}
                elif text == "END IONS":
                    if "TITLE" in current:
                        index[current["TITLE"]] = current
                elif "=" in text and text[:1].isalpha():
                    key, value = text.split("=", 1)
                    if key in ("TITLE", "PEPMASS", "CHARGE", "RTINSECONDS"):
                        current[key] = value
    return index


def validate_against_mgf(records: list[PrecursorRecord], mgf_dir: str) -> dict[str, object]:
    index = read_mgf_index(mgf_dir)
    matched = mz_mismatch = charge_mismatch = rt_mismatch = 0
    unmatched_examples = []
    for r in records:
        entry = index.get(r.spectrum_title)
        if entry is None:
            if len(unmatched_examples) < 5:
                unmatched_examples.append(r.spectrum_title)
            continue
        matched += 1
        pepmass = _num(entry.get("PEPMASS", "").split()[0] if entry.get("PEPMASS") else "")
        if pepmass is not None and r.precursor_mz is not None and abs(pepmass - r.precursor_mz) > 1e-4:
            mz_mismatch += 1
        charge_text = entry.get("CHARGE", "").rstrip("+")
        if charge_text.isdigit() and r.precursor_charge is not None and int(charge_text) != r.precursor_charge:
            charge_mismatch += 1
        rt_sec = _num(entry.get("RTINSECONDS", ""))
        if rt_sec is not None and r.retention_time is not None and abs(rt_sec / 60.0 - r.retention_time) > 1e-3:
            rt_mismatch += 1
    return {
        "mgf_dir": mgf_dir,
        "mgf_spectra": len(index),
        "inventory_spectra": len(records),
        "titles_matched": matched,
        "titles_unmatched": len(records) - matched,
        "mgf_titles_not_in_inventory": len(set(index) - {r.spectrum_title for r in records}),
        "precursor_mz_mismatches": mz_mismatch,
        "charge_mismatches": charge_mismatch,
        "retention_time_mismatches": rt_mismatch,
        "unmatched_examples": unmatched_examples,
    }


# --------------------------------------------------------------------------
# Orchestration
# --------------------------------------------------------------------------


def setup_logging(log_dir: str, verbosity: str) -> None:
    os.makedirs(log_dir, exist_ok=True)
    fmt_ = "%(asctime)s level=%(levelname)s logger=%(name)s %(message)s"
    root = logging.getLogger()
    root.handlers.clear()
    root.setLevel(logging.DEBUG)
    console = logging.StreamHandler(sys.stderr)
    console.setLevel(getattr(logging, verbosity.upper()))
    console.setFormatter(logging.Formatter(fmt_))
    root.addHandler(console)
    pipeline_file = logging.FileHandler(os.path.join(log_dir, "pipeline.log"), mode="w", encoding="utf-8")
    pipeline_file.setLevel(logging.DEBUG)
    pipeline_file.setFormatter(logging.Formatter(fmt_))
    root.addHandler(pipeline_file)
    extraction_file = logging.FileHandler(os.path.join(log_dir, "extraction.log"), mode="w", encoding="utf-8")
    extraction_file.setFormatter(logging.Formatter(fmt_))
    logging.getLogger("ms1quant.extraction").addHandler(extraction_file)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Identification-free MS1 precursor quantification of HCD MS2 spectra with FlashLFQ (Docker).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--project", default=None,
                   help="Project folder (e.g. /data/PXD009327). Derives --input <project>/MZML, "
                        "--output <project>/results and --reference-mgf-dir <project>/MGF (when it holds MGF files). "
                        "Explicit --input/--output/--reference-mgf-dir override the derived paths.")
    p.add_argument("--input", default=None, help="Directory with .mzML / .mzML.gz files (read-only).")
    p.add_argument("--output", default=None, help="Output directory (created if missing).")
    p.add_argument("--threads", type=int, default=4, help="CPU threads for FlashLFQ, Docker --cpus and extraction.")
    p.add_argument("--ppm", type=float, default=10.0, help="FlashLFQ peak-finding mass tolerance (ppm).")
    p.add_argument("--isotope-ppm", type=float, default=5.0, help="FlashLFQ isotope tolerance (ppm).")
    p.add_argument("--num-isotopes", type=int, default=2, help="FlashLFQ: isotope peaks required (>= 2).")
    p.add_argument("--id-specific-charge", action=argparse.BooleanOptionalAction, default=True,
                   help="FlashLFQ --chg: only search the precursor's own charge state. "
                        "--no-id-specific-charge searches all charge states (FlashLFQ's own default).")
    p.add_argument("--integrate", action=argparse.BooleanOptionalAction, default=True,
                   help="FlashLFQ --int: report envelope intensities summed over the peak. "
                        "--no-integrate reports the apex envelope only (FlashLFQ's own default).")
    p.add_argument("--image", default=DEFAULT_IMAGE, help="FlashLFQ Docker image (pinned by default).")
    p.add_argument("--docker-memory", default=None, help="Optional Docker --memory limit, e.g. 32g.")
    p.add_argument("--docker-arg", action="append", default=[], help="Extra argument passed to 'docker run' (repeatable).")
    p.add_argument("--timeout", type=float, default=None, help="Abort FlashLFQ after this many seconds.")
    p.add_argument("--tmp-dir", default=None, help="Directory for decompressed mzML staging (default: system temp).")
    p.add_argument("--keep-intermediate", action="store_true", help="Keep decompressed mzML staging files.")
    p.add_argument("--title-template", default="{run_id}: {native_id}",
                   help="Title used only when the mzML has no MS:1000796 spectrum title. "
                        "Placeholders: {run_id} {native_id} {file_stem} {scan_number} {index}.")
    p.add_argument("--reference-mgf-dir", default=None,
                   help="Optional directory with MGF files whose TITLEs should match the spectrum titles (validation only).")
    p.add_argument("--no-reference-mgf", action="store_true",
                   help="Do not use <project>/MGF for title validation even if it exists.")
    p.add_argument("--extract-only", action="store_true", help="Stop after writing the inventory and FlashLFQ input.")
    p.add_argument("--dry-run", action="store_true", help="Print the execution plan and exit.")
    p.add_argument("--log-level", default="info", choices=["debug", "info", "warning", "error"])
    # Development only: run a local FlashLFQ build instead of Docker.
    p.add_argument("--flashlfq-local-cmd", default=None, help=argparse.SUPPRESS)
    args = p.parse_args(argv)
    if args.threads < 1:
        p.error("--threads must be >= 1")
    if args.ppm <= 0 or args.isotope_ppm <= 0:
        p.error("tolerances must be > 0")
    if args.num_isotopes < 2:
        p.error("--num-isotopes must be >= 2 (FlashLFQ requirement)")
    try:
        resolve_project_paths(args)
    except PipelineError as exc:
        p.error(str(exc))
    return args


PROJECT_MZML_DIR = "MZML"
PROJECT_MGF_DIR = "MGF"
PROJECT_OUTPUT_DIR = "results"


def _find_subdir(parent: str, name: str) -> str | None:
    """Subdirectory of ``parent`` named ``name``; exact match first, then case-insensitive."""
    exact = os.path.join(parent, name)
    if os.path.isdir(exact):
        return exact
    for entry in sorted(os.listdir(parent)):
        if entry.lower() == name.lower() and os.path.isdir(os.path.join(parent, entry)):
            return os.path.join(parent, entry)
    return None


def resolve_project_paths(args: argparse.Namespace) -> None:
    """Fill --input/--output/--reference-mgf-dir from --project (PXD folder layout:
    <project>/MZML, <project>/MGF, results in <project>/results). Explicit options win."""
    if args.project:
        project = os.path.abspath(args.project)
        if not os.path.isdir(project):
            raise PipelineError(f"--project directory does not exist: {project}")
        if args.input is None:
            args.input = _find_subdir(project, PROJECT_MZML_DIR)
            if args.input is None:
                raise PipelineError(f"no {PROJECT_MZML_DIR}/ folder in {project}")
        if args.output is None:
            args.output = os.path.join(project, PROJECT_OUTPUT_DIR)
        if args.reference_mgf_dir is None and not args.no_reference_mgf:
            mgf_dir = _find_subdir(project, PROJECT_MGF_DIR)
            if mgf_dir and any(n.lower().endswith((".mgf", ".mgf.gz")) for n in os.listdir(mgf_dir)):
                args.reference_mgf_dir = mgf_dir
    if args.no_reference_mgf:
        args.reference_mgf_dir = None
    if args.input is None or args.output is None:
        raise PipelineError("give --project, or both --input and --output")


def summarize_values(values: list[float]) -> dict[str, float | int | None]:
    if not values:
        return {"n": 0}
    ordered = sorted(values)
    q = lambda f: ordered[min(len(ordered) - 1, int(f * (len(ordered) - 1)))]  # noqa: E731
    return {
        "n": len(ordered),
        "min": ordered[0],
        "q25": q(0.25),
        "median": statistics.median(ordered),
        "q75": q(0.75),
        "max": ordered[-1],
    }


def run_pipeline(args: argparse.Namespace) -> int:
    started = time.time()
    input_dir = os.path.abspath(args.input)
    output_dir = os.path.abspath(args.output)
    if output_dir == input_dir:
        raise PipelineError("--output must differ from --input")

    paths = discover_mzml_files(input_dir)
    params = FlashLfqParams(
        ppm=args.ppm,
        isotope_ppm=args.isotope_ppm,
        num_isotopes_required=args.num_isotopes,
        id_specific_charge=args.id_specific_charge,
        integrate=args.integrate,
        threads=args.threads,
    )
    planned_names = {os.path.basename(p): strip_mzml_suffix(os.path.basename(p)) + ".mzML" for p in paths}
    # docker rejects --cpus above the host CPU count
    docker_cpus = min(args.threads, os.cpu_count() or args.threads)

    if args.dry_run:
        staging_hint = os.path.join(os.path.abspath(args.tmp_dir) if args.tmp_dir else tempfile.gettempdir(),
                                    "ms1quant_staging_<random>")
        print("Execution plan (dry run; nothing is read or written)")
        print(f"  input directory : {input_dir}")
        print(f"  output directory: {output_dir}")
        print(f"  reference MGF   : {args.reference_mgf_dir or '(none)'}")
        print(f"  docker image    : {args.image}")
        print("  steps           : stream mzML -> precursor_inventory.tsv -> flashlfq_input.tsv -> one FlashLFQ"
              " container per run -> map QuantifiedPeaks.tsv -> ms1_precursor_quantification.tsv")
        for run_index, p in enumerate(paths, start=1):
            source = os.path.basename(p)
            run = strip_mzml_suffix(source)
            gz = is_gzip_file(p)
            cmd = build_docker_command(
                args.image, os.path.join(staging_hint, f"run{run_index:03d}"), planned_names[source], None if gz else p,
                os.path.join(output_dir, "flashlfq", run), params, docker_cpus, args.docker_memory, args.docker_arg,
            )
            print(f"  run {run}: {'decompress to staging, ' if gz else 'mount read-only, '}then")
            print("    " + shlex.join(cmd))
        return 0

    os.makedirs(output_dir, exist_ok=True)
    log_dir = os.path.join(output_dir, "logs")
    setup_logging(log_dir, args.log_level)
    xlog = logging.getLogger("ms1quant.extraction")
    log.info("event=start version=%s project=%s input=%s output=%s reference_mgf=%s files=%d", PIPELINE_VERSION,
             args.project, input_dir, output_dir, args.reference_mgf_dir, len(paths))

    checks: dict[str, dict[str, object]] = {}

    def check(name: str, passed: bool, detail: object = None) -> None:
        checks[name] = {"passed": bool(passed), "detail": detail}
        (log.info if passed else log.error)("event=check name=%s passed=%s detail=%s", name, passed, detail)

    # ---- 1. extraction ------------------------------------------------------
    t0 = time.time()
    extractions = extract_all(paths, args.title_template, args.threads)
    records: list[PrecursorRecord] = [r for e in extractions for r in e.records]
    for e in extractions:
        xlog.info(
            "event=extracted file=%s run_id=%s counts=%s title_sources=%s",
            e.source_file, e.run_id, json.dumps(e.spectrum_counts, sort_keys=True), json.dumps(e.title_source_counts),
        )
        status_counts = Counter(r.precursor_status for r in e.records)
        xlog.info("event=precursor_status file=%s counts=%s", e.source_file, json.dumps(dict(status_counts), sort_keys=True))
        mono_diff = sum(
            1 for r in e.records
            if r.thermo_monoisotopic_mz and r.precursor_mz and abs(r.thermo_monoisotopic_mz - r.precursor_mz) > 0.01
        )
        xlog.info("event=selected_vs_thermo_mono file=%s n_differ_gt_0.01=%d", e.source_file, mono_diff)
    log.info("event=extraction_done records=%d seconds=%.1f", len(records), time.time() - t0)

    check("hcd_spectra_found", len(records) > 0, {e.source_file: e.spectrum_counts.get("hcd_ms2", 0) for e in extractions})
    internal_ids = [r.internal_id for r in records]
    check("internal_ids_unique", len(set(internal_ids)) == len(internal_ids))
    titles = [r.spectrum_title for r in records]
    dup_titles = [t for t, n in Counter(titles).items() if n > 1]
    check("spectrum_titles_unique", not dup_titles, dup_titles[:5])
    check(
        "inventory_matches_spectrum_counts",
        all(e.spectrum_counts.get("hcd_ms2", 0) == len(e.records) for e in extractions),
    )

    assign_synthetic_ids(records)
    quantifiable = [r for r in records if r.precursor_status == PRECURSOR_OK]
    synthetic_ids = [r.synthetic_record_id for r in quantifiable]
    check("synthetic_ids_unique", len(set(synthetic_ids)) == len(synthetic_ids))

    inventory_path = os.path.join(output_dir, "precursor_inventory.tsv")
    write_tsv(inventory_path, INVENTORY_COLUMNS, ([getattr(r, c) for c in INVENTORY_COLUMNS] for r in records))
    log.info("event=written path=%s rows=%d", inventory_path, len(records))

    if args.reference_mgf_dir:
        mgf_report = validate_against_mgf(records, os.path.abspath(args.reference_mgf_dir))
        check(
            "titles_match_reference_mgf",
            mgf_report["titles_unmatched"] == 0
            and mgf_report["precursor_mz_mismatches"] == 0
            and mgf_report["charge_mismatches"] == 0
            and mgf_report["retention_time_mismatches"] == 0,
            mgf_report,
        )

    # ---- 2. FlashLFQ input ----------------------------------------------------
    staged_names = planned_names  # source_file -> FlashLFQ-visible file name
    input_path = os.path.join(output_dir, "flashlfq_input.tsv")
    n_sent = write_tsv(input_path, FLASHLFQ_INPUT_HEADER, flashlfq_input_rows(records, staged_names))
    per_file_sent = validate_flashlfq_input(input_path, set(staged_names.values()))
    check("flashlfq_input_valid_for_generic_reader", True, per_file_sent)
    check("flashlfq_input_rows_match_quantifiable", n_sent == len(quantifiable),
          {"sent": n_sent, "quantifiable": len(quantifiable)})
    log.info("event=written path=%s rows=%d", input_path, n_sent)

    summary: dict[str, object] = {
        "pipeline": {"name": "quantify_precursors.py", "version": PIPELINE_VERSION},
        "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(started)),
        "host": {"python": platform.python_version(), "platform": platform.platform()},
        "project_dir": os.path.abspath(args.project) if args.project else None,
        "input_dir": input_dir,
        "reference_mgf_dir": os.path.abspath(args.reference_mgf_dir) if args.reference_mgf_dir else None,
        "output_dir": output_dir,
        "parameters": {
            "ppm": args.ppm, "isotope_ppm": args.isotope_ppm, "num_isotopes_required": args.num_isotopes,
            "id_specific_charge_state": args.id_specific_charge, "integrate": args.integrate,
            "normalize": False, "match_between_runs": "neutralised (--pipfdr 0, one FlashLFQ run per file)",
            "threads": args.threads, "image": args.image, "title_template": args.title_template,
            "proton_mass": PROTON_MASS,
        },
        "files": [
            {
                "source_file": e.source_file,
                "run_id": e.run_id,
                "record_id_prefix": file_hash(e.source_file) + "_",
                "flashlfq_file_name": staged_names[e.source_file],
                "size_bytes": os.path.getsize(os.path.join(input_dir, e.source_file)),
                "spectrum_counts": e.spectrum_counts,
                "title_sources": e.title_source_counts,
                "precursor_status_counts": dict(Counter(r.precursor_status for r in e.records)),
            }
            for e in extractions
        ],
    }

    summary_path = os.path.join(output_dir, "quantification_summary.json")

    def write_summary() -> None:
        summary["checks"] = checks
        summary["all_checks_passed"] = all(c["passed"] for c in checks.values())
        summary["elapsed_seconds"] = round(time.time() - started, 1)
        with open(summary_path + ".part", "w", encoding="utf-8") as handle:
            json.dump(summary, handle, indent=2, sort_keys=False)
        os.replace(summary_path + ".part", summary_path)

    if args.extract_only:
        summary["mode"] = "extract_only"
        write_summary()
        log.info("event=done mode=extract_only")
        return 0 if summary["all_checks_passed"] else 2
    if n_sent == 0:
        raise PipelineError("No quantifiable HCD precursors; nothing to send to FlashLFQ")

    # ---- 3. FlashLFQ, one invocation per run ---------------------------------
    local_cmd = shlex.split(args.flashlfq_local_cmd) if args.flashlfq_local_cmd else None
    if local_cmd:
        log.warning("event=runner_local note=development mode, FlashLFQ is NOT running in Docker")
    else:
        check_docker(args.image)
    staging_parent = os.path.abspath(args.tmp_dir) if args.tmp_dir else None
    if staging_parent:
        os.makedirs(staging_parent, exist_ok=True)
    staging_root = tempfile.mkdtemp(prefix="ms1quant_staging_", dir=staging_parent)
    flash_root = os.path.join(output_dir, "flashlfq")
    records_by_file: dict[str, list[PrecursorRecord]] = defaultdict(list)
    for r in quantifiable:
        records_by_file[r.source_file].append(r)
    features: list[Feature] = []
    failed_files: dict[str, str] = {}
    run_reports: dict[str, dict[str, object]] = {}
    summary["flashlfq_runs"] = run_reports
    try:
        for run_index, path in enumerate(paths, start=1):
            source = os.path.basename(path)
            run = strip_mzml_suffix(source)
            run_records = records_by_file[source]
            report: dict[str, object] = {"runner": "local" if local_cmd else "docker", "records_sent": len(run_records)}
            run_reports[source] = report
            if not run_records:
                report["skipped"] = "no quantifiable precursors"
                continue
            work_dir = os.path.join(flash_root, run)
            out_dir = os.path.join(work_dir, RUN_OUTPUT_DIR)
            os.makedirs(out_dir, exist_ok=True)
            for name in FLASHLFQ_OUTPUT_FILES:  # stale results from an earlier run
                if os.path.isfile(os.path.join(out_dir, name)):
                    os.remove(os.path.join(out_dir, name))
            run_input = os.path.join(work_dir, RUN_INPUT_NAME)
            write_tsv(run_input, FLASHLFQ_INPUT_HEADER, flashlfq_input_rows(run_records, staged_names))
            spectra_dir = os.path.join(staging_root, f"run{run_index:03d}")  # short: path length limits
            try:
                name, readonly_source = stage_spectra_file(path, spectra_dir, copy_plain=bool(local_cmd))
                if name != staged_names[source]:
                    raise PipelineError(f"staged name {name} differs from {staged_names[source]}")
                kill_cmd = None
                if local_cmd:
                    cmd = local_cmd + params.cli_args(idt=run_input, rep=spectra_dir, out=out_dir)
                else:
                    container = f"ms1quant-{os.getpid()}-{run_index:03d}"
                    cmd = build_docker_command(
                        args.image, spectra_dir, name, readonly_source, work_dir, params,
                        docker_cpus, args.docker_memory, args.docker_arg, container_name=container,
                    )
                    kill_cmd = ["docker", "kill", container]
                t1 = time.time()
                rc, console = run_logged(cmd, os.path.join(log_dir, f"flashlfq_{run}.log"), args.timeout, kill_cmd)
            except (PipelineError, OSError) as exc:
                failed_files[source] = str(exc)
                report["error"] = str(exc)
                log.error("event=flashlfq_failed file=%s error=%s", source, exc)
                continue
            finally:
                if not args.keep_intermediate:
                    shutil.rmtree(spectra_dir, ignore_errors=True)
            console_report = check_flashlfq_console(console, len(run_records))
            settings_path = os.path.join(out_dir, "FlashLfqSettings.toml")
            settings = read_flashlfq_settings(settings_path) if os.path.isfile(settings_path) else {}
            mismatch = (settings_mismatches(settings, params.expected_settings())
                        if settings else {"FlashLfqSettings.toml": "missing"})
            peaks_path = os.path.join(out_dir, "QuantifiedPeaks.tsv")
            report.update(
                command=cmd, exit_code=rc, seconds=round(time.time() - t1, 1),
                console_checks=console_report, effective_settings=settings, settings_mismatches=mismatch,
            )
            log.info("event=flashlfq_done file=%s rc=%d report=%s", source, rc, json.dumps(console_report))
            if rc != 0 or console_report["failure_markers"] or mismatch or not os.path.isfile(peaks_path):
                reason = (f"exit code {rc}; failure markers {console_report['failure_markers']}; "
                          f"settings mismatches {mismatch}; QuantifiedPeaks.tsv present: {os.path.isfile(peaks_path)}")
                failed_files[source] = reason
                report["error"] = reason
                log.error("event=flashlfq_failed file=%s reason=%s log=%s", source, reason,
                          os.path.join(log_dir, f"flashlfq_{run}.log"))
                continue
            features.extend(parse_quantified_peaks(peaks_path))
    finally:
        if args.keep_intermediate:
            log.info("event=kept_staging path=%s", staging_root)
        else:
            shutil.rmtree(staging_root, ignore_errors=True)

    check("flashlfq_runs_succeeded", not failed_files, failed_files)
    ok_reports = {k: v for k, v in run_reports.items() if "console_checks" in v and k not in failed_files}
    check(
        "flashlfq_read_all_records",
        all(v["console_checks"]["records_reported_read"] == v["records_sent"]
            and v["console_checks"]["reader_rejection_messages"] == 0 for v in ok_reports.values()),
        {k: v["console_checks"] for k, v in ok_reports.items()},
    )
    check(
        "flashlfq_effective_settings_as_requested",
        all(not v.get("settings_mismatches") for v in run_reports.values() if "settings_mismatches" in v),
        {k: v.get("settings_mismatches") for k, v in run_reports.items()},
    )
    mbr_rows = sum(1 for f in features if "MBR" in f.detection_type.upper())
    check("no_match_between_runs_peaks", mbr_rows == 0,
          {"mbr_rows": mbr_rows, "detection_types": dict(Counter(f.detection_type for f in features))})

    # ---- 4. mapping -----------------------------------------------------------
    assign_feature_ids(features)
    mapping, problems = map_features(records, features, staged_names)
    check("feature_mapping_integrity", not problems, problems[:20])
    check("no_cross_file_assignments", not any(p.startswith("cross-file") for p in problems))
    missing = [r.synthetic_record_id for r in quantifiable
               if r.source_file not in failed_files and r.synthetic_record_id not in mapping]
    check("every_record_in_flashlfq_output", not missing, {"missing": len(missing), "examples": missing[:5]})

    # round-trip: single-record peaks echo the record's mass/charge/RT
    record_by_id = {r.synthetic_record_id: r for r in quantifiable}
    roundtrip_bad = 0
    for feature in features:
        if len(feature.record_ids) == 1 and feature.record_ids[0] in record_by_id:
            r = record_by_id[feature.record_ids[0]]
            if (
                feature.first_record_charge != r.precursor_charge
                or feature.first_record_mass is None
                or abs(feature.first_record_mass - r.neutral_mass) > 1e-5
                or (feature.first_record_rt is not None and abs(feature.first_record_rt - r.retention_time) > 1e-5)
            ):
                roundtrip_bad += 1
    check("flashlfq_echoes_input_values", roundtrip_bad == 0, {"mismatching_single_record_peaks": roundtrip_bad})

    # ---- 5. outputs -------------------------------------------------------------
    final_rows = []
    mapping_rows = []
    for r in records:
        feature = mapping.get(r.synthetic_record_id) if r.synthetic_record_id else None
        if r.precursor_status != PRECURSOR_OK:
            status = f"not_quantified_{r.precursor_status}"
        elif r.source_file in failed_files:
            status = STATUS_RUN_FAILED
        else:
            status = quant_status(feature)
        row = final_row(r, feature, status)
        final_rows.append([row.get(c) for c in FINAL_COLUMNS])
        if feature is not None:
            mapping_rows.append([
                r.internal_id, r.source_file, r.spectrum_title, r.native_scan_id, r.synthetic_record_id,
                feature.feature_id or None, feature.file_name,
                len(feature.record_ids) if has_apex(feature) else None,
                feature.rt_apex, feature.peak_mz, feature.intensity if has_apex(feature) else None,
                status, "record_id_in_QuantifiedPeaks_Full_Sequence",
            ])

    final_path = os.path.join(output_dir, "ms1_precursor_quantification.tsv")
    n_final = write_tsv(final_path, FINAL_COLUMNS, final_rows)
    mapping_path = os.path.join(output_dir, "feature_mapping.tsv")
    n_mapping = write_tsv(mapping_path, MAPPING_COLUMNS, mapping_rows)
    log.info("event=written path=%s rows=%d", final_path, n_final)
    log.info("event=written path=%s rows=%d", mapping_path, n_mapping)

    check("one_output_row_per_hcd_spectrum", n_final == len(records), {"rows": n_final, "hcd_spectra": len(records)})
    status_idx = FINAL_COLUMNS.index("quantification_status")
    per_file_out = Counter(row[0] for row in final_rows)
    check(
        "per_file_output_counts_match_inventory",
        all(per_file_out[e.source_file] == len(e.records) for e in extractions),
        dict(per_file_out),
    )

    quantified_features = [f for f in features if has_apex(f)]
    members = Counter(len(f.record_ids) for f in quantified_features)
    statuses = Counter(row[status_idx] for row in final_rows)
    per_file_status: dict[str, Counter[str]] = defaultdict(Counter)
    for row in final_rows:
        per_file_status[row[0]][row[status_idx]] += 1
    charge_match = Counter(
        row[FINAL_COLUMNS.index("feature_charge_matches_precursor")] for row in final_rows
        if row[status_idx] == STATUS_QUANTIFIED
    )
    rt_within = Counter(
        row[FINAL_COLUMNS.index("ms2_rt_within_feature_bounds")] for row in final_rows
        if row[status_idx] == STATUS_QUANTIFIED
    )
    summary["counts"] = {
        "hcd_ms2_spectra": len(records),
        "quantifiable_precursors": len(quantifiable),
        "records_sent_to_flashlfq": n_sent,
        "flashlfq_failed_files": sorted(failed_files),
        "flashlfq_peak_rows": len(features),
        "features_with_apex": len(quantified_features),
        "features_shared_by_multiple_spectra": sum(n for k, n in members.items() if k > 1),
        "spectra_in_shared_features": sum(k * n for k, n in members.items() if k > 1),
        "max_spectra_per_feature": max(members) if members else 0,
        "feature_size_distribution": {str(k): members[k] for k in sorted(members)},
        "quantification_status": dict(statuses),
        "quantification_status_per_file": {k: dict(v) for k, v in per_file_status.items()},
        "feature_charge_matches_precursor_charge": {fmt(k): v for k, v in charge_match.items()},
        "ms2_rt_within_feature_bounds": {fmt(k): v for k, v in rt_within.items()},
        "abs_ms2_rt_minus_feature_apex_min": summarize_values(
            [abs(row[FINAL_COLUMNS.index("ms2_rt_minus_feature_apex")]) for row in final_rows
             if row[status_idx] == STATUS_QUANTIFIED]
        ),
        "ms1_intensity": summarize_values(
            [row[FINAL_COLUMNS.index("ms1_intensity")] for row in final_rows if row[status_idx] == STATUS_QUANTIFIED]
        ),
    }
    write_summary()
    log.info("event=done status_counts=%s all_checks_passed=%s", json.dumps(dict(statuses)), summary["all_checks_passed"])
    return 0 if summary["all_checks_passed"] else 2


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        return run_pipeline(args)
    except PipelineError as exc:
        if logging.getLogger().handlers:
            log.error("event=fatal error=%s", exc)
        else:
            print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
