"""Stand-in for the FlashLFQ CLI, used to test orchestration without Docker.

It reads the Generic identification file like FlashLFQ, writes a
FlashLfqSettings.toml reflecting the switches it was given, and writes a
QuantifiedPeaks.tsv with the same columns as FlashLFQ 2.1.5:
- the first two records are merged into one peak ("A|B", PSMs Mapped 2),
- the third record has no apex (intensity 0, "-" fields),
- every other record gets its own peak.
Set FAKE_FLASHLFQ_FAIL=1 to simulate FlashLFQ's silent failure (exit 0) and
FAKE_FLASHLFQ_SLEEP=<seconds> to simulate a hang.
"""

import csv
import os
import sys
import time

HEADER = [
    "File Name", "Base Sequence", "Full Sequence", "Protein Group", "Organism",
    "Peptide Monoisotopic Mass", "MS2 Retention Time", "Precursor Charge", "Theoretical MZ",
    "Peak intensity", "Peak RT Start", "Peak RT Apex", "Peak RT End", "Peak FWHM", "Peak FWHM Status",
    "Peak MZ", "Peak Charge", "Num Charge States Observed", "Peak Detection Type", "PIP Q-Value",
    "PIP PEP", "PSMs Mapped", "Base Sequences Mapped", "Full Sequences Mapped",
    "Peak Split Valley RT", "Peak Apex Mass Error (ppm)", "Decoy Peptide", "Random RT",
]


def main(argv):
    args = {}
    i = 0
    while i < len(argv):
        token = argv[i]
        if i + 1 < len(argv) and not argv[i + 1].startswith("--"):
            args[token] = argv[i + 1]
            i += 2
        else:
            args[token] = True
            i += 1
    if os.environ.get("FAKE_FLASHLFQ_SLEEP"):
        time.sleep(float(os.environ["FAKE_FLASHLFQ_SLEEP"]))
    if os.environ.get("FAKE_FLASHLFQ_FAIL"):
        print("Problem reading PSMs: simulated")
        return 0
    with open(args["--idt"], newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    print(f"Done reading PSMs; found {len(rows)}")
    print(f"Setup is OK; read in {len(rows)} identifications; starting FlashLFQ engine")
    out = args["--out"]
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "FlashLfqSettings.toml"), "w") as toml:
        toml.write(
            "Silent = false\nNormalize = {nor}\nPpmTolerance = {ppm}\nIsotopePpmTolerance = {iso}\n"
            "Integrate = {int_}\nNumIsotopesRequired = {nis}\nIdSpecificChargeState = {chg}\n"
            "MatchBetweenRuns = true\nBayesianProteinQuant = false\nUseSharedPeptidesForProteinQuant = false\n"
            "MbrQValueThreshold = {pip}\nUsePepQValue = false\nRnaMode = false\n".format(
                nor="true" if "--nor" in args else "false",
                ppm=float(args["--ppm"]), iso=float(args["--iso"]),
                int_="true" if "--int" in args else "false", nis=args["--nis"],
                chg="true" if "--chg" in args else "false", pip=float(args.get("--pipfdr", 0.05)),
            )
        )
    peaks = []
    if len(rows) >= 2:
        peaks.append((rows[:2], True))
        rest = rows[2:]
    else:
        rest = rows
    if rest:
        peaks.append(([rest[0]], False))
        peaks.extend(([r], True) for r in rest[1:])
    with open(os.path.join(out, "QuantifiedPeaks.tsv"), "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(HEADER)
        for members, apex in peaks:
            first = members[0]
            rt = float(first["Scan Retention Time"])
            z = int(first["Precursor Charge"])
            mass = float(first["Peptide Monoisotopic Mass"])
            mz = mass / z + 1.007276466879
            file_name = os.path.splitext(first["File Name"])[0]
            ids = "|".join(m["Full Sequence"] for m in members)
            if apex:
                fields = ["1000.5", repr(rt - 0.1), repr(rt), repr(rt + 0.1), "0.05", "Measured", repr(mz), str(z), "1"]
            else:
                fields = ["0", "-", "-", "-", "-", "NoApex", "-", "-", "0"]
            writer.writerow([
                file_name, ids, ids, "NO_PROTEIN_INFERENCE", "", first["Peptide Monoisotopic Mass"],
                first["Scan Retention Time"], first["Precursor Charge"], repr(mz),
                *fields[:7], fields[7], fields[8], "MSMS", "", "", str(len(members)), str(len(members)),
                str(len(members)), "0", "0.5" if apex else "NaN", "False", "False",
            ])
    print("Done quantifying")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
