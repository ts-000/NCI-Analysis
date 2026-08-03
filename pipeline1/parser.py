import csv
import re
import sys
from pathlib import Path

TABLE = "second order perturbation theory analysis of fock matrix in nbo basis"
END = "natural bond orbitals (summary):"
ROW = re.compile(
    r"^\s*(\d+\.)\s*(.+?)\s+(\d+\.)\s*(.+?)\s+"
    r"([0-9]*\.[0-9]+)\s+([0-9]*\.[0-9]+)\s+([0-9]*\.[0-9]+)\s*$"
)


def parse_nbo_e2(text):
    lines = text.splitlines()
    end = next(
        (i for i, line in enumerate(lines) if "Normal termination of Gaussian" in line),
        len(lines),
    )
    # Use the last E(2) section in the main Gaussian job. This is the section corresponding to the last geometry selected by molecule_visualizer
    tables = [i for i, line in enumerate(lines[:end]) if TABLE in line.lower()]
    if not tables:
        return []
    rows = []
    for line in lines[tables[-1] + 1:end]:
        if END in line.lower():
            break
        match = ROW.match(line)
        if match:
            donor_number, donor, acceptor_number, acceptor, e2, _, _ = match.groups()
            rows.append((
                f"{donor_number} {donor.strip()}",
                f"{acceptor_number} {acceptor.strip()}",
                e2,
            ))
    return rows


def write_csv(rows, outname):
    with open(outname, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["donor", "acceptor", "E2_kcal/mol"])
        writer.writerows(rows)


def main():
    if len(sys.argv) < 2:
        print("Usage: python parser.py file1.log [file2.log ...]")
        return 1
    failed = False
    for infile in sys.argv[1:]:
        try:
            text = Path(infile).read_text(encoding="utf-8", errors="ignore")
            rows = parse_nbo_e2(text)
        except Exception as exc:
            print(f"Error reading {infile}: {exc}")
            failed = True
            continue
        if not rows:
            print("No E(2) entries found in", infile)
            failed = True
            continue
        outcsv = str(Path(infile).with_suffix(".csv"))
        write_csv(rows, outcsv)
        print("Wrote", len(rows), "rows to", outcsv)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
