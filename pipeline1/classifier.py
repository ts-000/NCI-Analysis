import csv
import re
import sys
from pathlib import Path


def extract_orbital(label):
    match = re.search(r"(LP|BD\*?)\s*\(\s*(\d+)\s*\)", label.upper())
    return (match.group(1), int(match.group(2))) if match else (None, None)


def classify_interaction(donor, acceptor):
    donor_type, donor_number = extract_orbital(donor)
    acceptor_type, acceptor_number = extract_orbital(acceptor)
    if donor_type == "LP" and acceptor_type == "BD*":
        return "n→σ*" if acceptor_number == 1 else "n→π*" if acceptor_number == 2 else None
    if donor_type == "BD" and acceptor_type == "BD*":
        donor_name = "σ" if donor_number == 1 else "π" if donor_number == 2 else None
        acceptor_name = "σ*" if acceptor_number == 1 else "π*" if acceptor_number == 2 else None
        return f"{donor_name}→{acceptor_name}" if donor_name and acceptor_name else None
    return None


def main():
    if len(sys.argv) < 2:
        print("Usage: python classifier.py filtered_file1.csv [filtered_file2.csv ...]")
        return 1
    headers = [
        "donor", "acceptor", "E2_kcal/mol",
        "donor_labels", "acceptor_labels", "bond_distance",
        "interaction_type",
    ]
    for infile in sys.argv[1:]:
        inpath = Path(infile)
        base = inpath.name[len("filtered_"):] if inpath.name.startswith("filtered_") else inpath.name
        outfile = inpath.parent / f"classified_{base}"
        with inpath.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        out_rows = []
        for row in rows:
            interaction_type = classify_interaction(row["donor"], row["acceptor"])
            if interaction_type is None:
                continue
            out_rows.append({
                "donor": row["donor"],
                "acceptor": row["acceptor"],
                "E2_kcal/mol": row["E2_kcal/mol"],
                "donor_labels": row["donor_labels"],
                "acceptor_labels": row["acceptor_labels"],
                "bond_distance": row["bond_distance"],
                "interaction_type": interaction_type,
            })
        with outfile.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=headers)
            writer.writeheader()
            writer.writerows(out_rows)
        print(f"Wrote {len(out_rows)} rows to {outfile}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
