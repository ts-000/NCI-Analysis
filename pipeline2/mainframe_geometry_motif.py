from pathlib import Path
import re
import numpy as np
import pandas as pd


OUTFILE = "mainframe_biomotifs.csv"
HERE = Path(__file__).resolve().parent

NEED = {
    "donor",
    "acceptor",
    "donor_labels",
    "acceptor_labels",
    "E2_kcal/mol",
    "bond_distance",
    "interaction_type",
}

OUT_COLS = [
    "interaction_id",
    "motif_id",
    "E2_kcal/mol",
    "donor_orbital",
    "acceptor_orbital",
    "donor_labels",
    "acceptor_labels",
    "bond_distance",
    "interaction_type",
    "distance",
    "angle",
    "dihedral",
]


def dist(a, b):
    return float(np.linalg.norm(a - b))


def vec_angle(u, v):
    nu, nv = np.linalg.norm(u), np.linalg.norm(v)
    if nu == 0 or nv == 0:
        return np.nan
    cs = np.clip(np.dot(u, v) / (nu * nv), -1.0, 1.0)
    return float(np.degrees(np.arccos(cs)))


def dihedral(p0, p1, p2, p3):
    b0 = -(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2
    b1n = np.linalg.norm(b1)
    if b1n == 0:
        return np.nan
    b1 = b1 / b1n
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    if np.linalg.norm(v) == 0 or np.linalg.norm(w) == 0:
        return np.nan
    return float(
        np.degrees(
            np.arctan2(
                np.dot(np.cross(b1, v), w),
                np.dot(v, w),
            )
        )
    )


def abs_dih(p0, p1, p2, p3):
    v = dihedral(p0, p1, p2, p3)
    return abs(v) if not pd.isna(v) else np.nan


def load_xyz(fp):
    with open(fp, encoding="utf-8") as fh:
        lns = fh.readlines()
    natom = int(lns[0].strip())
    xyz = {}
    for i, line in enumerate(lns[2:2 + natom], start=1):
        atom, x, y, z = line.split()[:4]
        xyz[f"{atom}{i}"] = np.array(
            [float(x), float(y), float(z)], dtype=float
        )
    return xyz


def atom_list(v):
    return re.findall(r"[A-Z][a-z]?\d+", str(v))


def tidy_atoms(v):
    return "|".join(atom_list(v))


def no_serial(v):
    if pd.isna(v):
        return v
    return re.sub(r"^\s*\d+\.\s*", "", str(v))


def clean_orb(v):
    txt = re.sub(r"\s+", " ", str(no_serial(v))).strip()
    txt = re.sub(r"\(\s*(\d+)\s*\)", r"(\1)", txt)
    txt = re.sub(
        r"\b([A-Z][a-z]?)\s*(\d+)\b",
        lambda match: f"{match.group(1)}{int(match.group(2))}",
        txt,
    )
    return re.sub(r"\s*-\s*", "-", txt)


def norm_kind(v):
    txt = str(v).strip().lower()
    txt = txt.replace("σ", "sig").replace("sigma", "sig")
    txt = txt.replace("π", "pi")
    txt = txt.replace("→", "-").replace("->", "-")
    txt = txt.replace("–", "-").replace("—", "-")
    txt = txt.replace("∗", "").replace("*", "")
    txt = re.sub(r"[\s_]+", "", txt)
    return re.sub(r"-+", "-", txt)


def elem(lbl):
    match = re.match(r"[A-Z][a-z]?", str(lbl))
    if not match:
        raise ValueError(f"Cannot determine element from atom label: {lbl}")
    return match.group(0)


def elems_are(labels, want):
    return sorted(elem(lbl) for lbl in labels) == sorted(want)


def get_atom(labels, element):
    for lbl in labels:
        if elem(lbl) == element:
            return lbl
    raise ValueError(f"Expected element {element} in labels: {'|'.join(labels)}")


def get_two(labels, element):
    got = [lbl for lbl in labels if elem(lbl) == element]
    if len(got) != 2:
        raise ValueError(f"Expected two {element} atoms in labels: {'|'.join(labels)}")
    return got


def atom_dist(label_a, label_b, xyz):
    return dist(xyz[label_a], xyz[label_b])


def near_far(ref, a1, a2, xyz):
    if atom_dist(ref, a1, xyz) <= atom_dist(ref, a2, xyz):
        return a1, a2
    return a2, a1


def donor_carbon(da, ac, xyz):
    donor_element = elem(da)
    lim = 1.65 if donor_element == "O" else 1.75
    cand = [
        atom
        for atom in xyz
        if atom != ac
        and elem(atom) == "C"
        and 1.10 < atom_dist(da, atom, xyz) < lim
    ]
    if not cand:
        raise ValueError(
            f"No carbon bonded directly to donor atom {da}"
        )
    if len(cand) == 1:
        return cand[0]
    return min(
        cand,
        key=lambda atom: atom_dist(atom, ac, xyz),
    )


def geom_out(dv=np.nan, av=np.nan, tv=np.nan):
    return {
        "distance": dv,
        "angle": av,
        "dihedral": tv,
    }


def geom_npi(donor, acceptor, xyz):
    da = donor[0]
    ac = get_atom(acceptor, "C")
    ao = get_atom(acceptor, "O")
    return geom_out(
        dv=atom_dist(da, ac, xyz),
        av=vec_angle(
            xyz[da] - xyz[ac],
            xyz[ao] - xyz[ac],
        ),
    )


def geom_nsig_ch(donor, acceptor, xyz):
    da = donor[0]
    ac = get_atom(acceptor, "C")
    ah = get_atom(acceptor, "H")
    dc = donor_carbon(
        da,
        ac,
        xyz,
    )
    return geom_out(
        dv=atom_dist(da, ah, xyz),
        tv=abs_dih(
            xyz[da],
            xyz[dc],
            xyz[ac],
            xyz[ah],
        ),
    )


def geom_nsig_cc(donor, acceptor, xyz):
    da = donor[0]
    c1, c2 = get_two(acceptor, "C")
    ca, cb = near_far(da, c1, c2, xyz)
    return geom_out(
        av=vec_angle(
            xyz[da] - xyz[ca],
            xyz[cb] - xyz[ca],
        )
    )


def geom_pipi(donor, acceptor, xyz):
    dc = get_atom(donor, "C")
    do = get_atom(donor, "O")
    ac = get_atom(acceptor, "C")
    ao = get_atom(acceptor, "O")
    return geom_out(
        dv=atom_dist(do, ao, xyz),
        tv=abs_dih(
            xyz[do],
            xyz[dc],
            xyz[ac],
            xyz[ao],
        ),
    )


def geom_sigpi_nc(donor, acceptor, xyz):
    dc = get_atom(donor, "C")
    dn = get_atom(donor, "N")
    ac = get_atom(acceptor, "C")
    ao = get_atom(acceptor, "O")
    return geom_out(
        dv=atom_dist(dc, ao, xyz),
        tv=dihedral(
            xyz[dc],
            xyz[dn],
            xyz[ac],
            xyz[ao],
        ),
    )


def geom_sigpi_chco(donor, acceptor, xyz):
    dc = get_atom(donor, "C")
    dh = get_atom(donor, "H")
    ac = get_atom(acceptor, "C")
    ao = get_atom(acceptor, "O")
    return geom_out(
        dv=atom_dist(dh, ao, xyz),
        tv=dihedral(
            xyz[dh],
            xyz[dc],
            xyz[ac],
            xyz[ao],
        ),
    )


def geom_sigpi_chcc(donor, acceptor, xyz):
    dc = get_atom(donor, "C")
    dh = get_atom(donor, "H")
    c1, c2 = get_two(acceptor, "C")
    hit_c, other_c = near_far(dh, c1, c2, xyz)
    return geom_out(
        dv=atom_dist(dh, hit_c, xyz),
        av=vec_angle(
            xyz[dc] - xyz[hit_c],
            xyz[other_c] - xyz[hit_c],
        ),
    )


def geom_pisig(donor, acceptor, xyz):
    do = get_atom(donor, "O")
    ac = get_atom(acceptor, "C")
    ah = get_atom(acceptor, "H")
    return geom_out(
        av=vec_angle(
            xyz[ah] - xyz[ac],
            xyz[do] - xyz[ac],
        )
    )


def geom_sigsig_cc(donor, acceptor, xyz):
    dc = get_atom(donor, "C")
    dh = get_atom(donor, "H")
    c1, c2 = get_two(acceptor, "C")
    near_c, far_c = near_far(
        dc, c1, c2, xyz
    )
    return geom_out(
        dv=atom_dist(dc, far_c, xyz),
        tv=dihedral(
            xyz[dh],
            xyz[dc],
            xyz[near_c],
            xyz[far_c],
        ),
    )


def geom_sigsig_ch(donor, acceptor, xyz):
    dc = get_atom(donor, "C")
    dh = get_atom(donor, "H")
    ac = get_atom(acceptor, "C")
    ah = get_atom(acceptor, "H")
    return geom_out(
        dv=atom_dist(dh, ah, xyz),
        tv=dihedral(
            xyz[dh],
            xyz[dc],
            xyz[ac],
            xyz[ah],
        ),
    )


def geom_sigsig_cn(donor, acceptor, xyz):
    dc = get_atom(donor, "C")
    dh = get_atom(donor, "H")
    ac = get_atom(acceptor, "C")
    an = get_atom(acceptor, "N")
    return geom_out(
        dv=atom_dist(dh, ac, xyz),
        tv=dihedral(
            xyz[dh],
            xyz[dc],
            xyz[an],
            xyz[ac],
        ),
    )


def calc_geom(row, xyz):
    itype = norm_kind(row["interaction_type"])
    donor = atom_list(row["donor_labels"])
    acceptor = atom_list(row["acceptor_labels"])
    blank = geom_out()
    missing_atoms = [atom for atom in donor + acceptor if atom not in xyz]
    if missing_atoms:
        raise KeyError(
            "Atoms not found in XYZ coordinates: " + ", ".join(missing_atoms)
        )
    if itype == "n-pi":
        if (
            len(donor) == 1
            and elem(donor[0]) in {"N", "O"}
            and elems_are(acceptor, ("C", "O"))
        ):
            return geom_npi(donor, acceptor, xyz)
        return blank

    if itype in {"n-sig", "n-sigch", "n-sigcc"}:
        valid_donor = (
            len(donor) == 1
            and elem(donor[0]) in {"N", "O"}
        )
        if not valid_donor:
            return blank

        if (
            itype != "n-sigcc"
            and elems_are(acceptor, ("C", "H"))
        ):
            return geom_nsig_ch(donor, acceptor, xyz)

        if (
            itype != "n-sigch"
            and elems_are(acceptor, ("C", "C"))
        ):
            return geom_nsig_cc(donor, acceptor, xyz)

        return blank

    if itype == "pi-pi":
        if elems_are(donor, ("C", "O")) and elems_are(
            acceptor, ("C", "O")
        ):
            return geom_pipi(donor, acceptor, xyz)
        return blank

    if itype == "sig-pi":
        if elems_are(donor, ("C", "N")) and elems_are(
            acceptor, ("C", "O")
        ):
            return geom_sigpi_nc(donor, acceptor, xyz)
        if elems_are(donor, ("C", "H")) and elems_are(
            acceptor, ("C", "O")
        ):
            return geom_sigpi_chco(donor, acceptor, xyz)
        if elems_are(donor, ("C", "H")) and elems_are(
            acceptor, ("C", "C")
        ):
            return geom_sigpi_chcc(donor, acceptor, xyz)
        return blank

    if itype == "pi-sig":
        if elems_are(donor, ("C", "O")) and elems_are(
            acceptor, ("C", "H")
        ):
            return geom_pisig(donor, acceptor, xyz)
        return blank

    if itype == "sig-sig":
        if not elems_are(donor, ("C", "H")):
            return blank
        if elems_are(acceptor, ("C", "C")):
            return geom_sigsig_cc(donor, acceptor, xyz)
        if elems_are(acceptor, ("C", "H")):
            return geom_sigsig_ch(donor, acceptor, xyz)
        if elems_are(acceptor, ("C", "N")):
            return geom_sigsig_cn(donor, acceptor, xyz)
        return blank

    return blank


def make_id(row):
    return "__".join(
        str(row[column])
        for column in (
            "motif_id",
            "donor_orbital",
            "acceptor_orbital",
            "donor_labels",
            "acceptor_labels",
            "interaction_type",
        )
    )


def load_motif(mid, csvf, xyzf):
    df = pd.read_csv(csvf)
    missing = NEED - set(df.columns)
    if missing:
        raise ValueError(
            f"{csvf.name} is missing required columns: "
            + ", ".join(sorted(missing))
        )

    xyz = load_xyz(xyzf)
    df["donor_orbital"] = df["donor"].apply(clean_orb)
    df["acceptor_orbital"] = df["acceptor"].apply(clean_orb)
    df["donor_labels"] = df["donor_labels"].apply(tidy_atoms)
    df["acceptor_labels"] = df["acceptor_labels"].apply(tidy_atoms)
    df["E2_kcal/mol"] = pd.to_numeric(df["E2_kcal/mol"])

    geom_rows = []
    for row_no, (_, row) in enumerate(df.iterrows(), start=2):
        try:
            geom_rows.append(calc_geom(row, xyz))
        except Exception as exc:
            raise ValueError(
                f"Geometry error in motif '{mid}', "
                f"CSV row {row_no}: {exc}"
            ) from exc

    if geom_rows:
        geom_df = pd.DataFrame(geom_rows, index=df.index)
        df[["distance", "angle", "dihedral"]] = geom_df[
            ["distance", "angle", "dihedral"]
        ]
    else:
        df[["distance", "angle", "dihedral"]] = np.nan

    df["motif_id"] = mid
    df["interaction_id"] = df.apply(make_id, axis=1)
    return df[OUT_COLS]


# pairs use whatever comes after classified_
def find_pairs():
    csv_files = sorted(HERE.glob("classified_*.csv"))
    if not csv_files:
        raise ValueError("No classified_*.csv motif files were found beside the script.")

    return [
        (
            csvf.stem.removeprefix("classified_"),
            csvf,
            HERE / f"{csvf.stem.removeprefix('classified_')}.xyz",
        )
        for csvf in csv_files
    ]


def main():
    chunks = []
    for mid, csvf, xyzf in find_pairs():
        print(f"Processing motif: {mid}")
        chunks.append(load_motif(mid, csvf, xyzf))

    out = pd.concat(chunks, ignore_index=True)
    if out.empty:
        raise ValueError("No active interaction rows were found.")

    out.to_csv(HERE / OUTFILE, index=False)
    print(f"DONE: {OUTFILE} created")
    print(f"Motifs processed: {out['motif_id'].nunique()}")
    print(f"Active interactions: {len(out)}")


if __name__ == "__main__":
    main()
