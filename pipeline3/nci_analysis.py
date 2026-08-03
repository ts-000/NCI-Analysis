import os
import re
import math
import warnings
import shutil
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from mpl_toolkits.mplot3d import proj3d

warnings.filterwarnings("ignore")

SCRIPT_VERSION = "nci_analysis_v9_h"
print(f"RUNNING {SCRIPT_VERSION}")

FILE = "mainframe.csv"
MOTIF_FILE = "mainframe_biomotifs.csv"
OUTDIR = "mega_nci_report"
FIGDIR = os.path.join(OUTDIR, "figures")
FIGPRESENTDIR = os.path.join(OUTDIR, "figures_presentation")
FIGPAPERDIR = os.path.join(OUTDIR, "figures_paper")
FIGSUPPDIR = os.path.join(OUTDIR, "figures_supplementary")
PAPERTABLEDIR = os.path.join(OUTDIR, "tables_paper")
SUPPTABLEDIR = os.path.join(OUTDIR, "tables_supplementary")

CONFORMER_ORDER = ["N1", "N2", "P1", "P2"]

HAMMETT_SIGMA_P = {
    "NMe2": -0.83,
    "NH2": -0.66,
    "OH": -0.37,
    "H": 0.00,
    "F": 0.06,
    "Cl": 0.23,
    "CF3": 0.54,
    "CN": 0.66,
    "NO2": 0.78,
}
SUB_ORDER = sorted(HAMMETT_SIGMA_P, key=HAMMETT_SIGMA_P.get)

MIN_CORR_N = 6
LINE_DETAIL_R_THRESHOLD = 0.75

CORE_MAX_ATOM = 27

EXPORT_EXCEL_PLAYGROUND = True
EXCEL_PLAYGROUND_NAME = "nci_energy_modular_table.xlsx"

PAPER_TABLE_WORKBOOK = "paper_figure_tables.xlsx"

shutil.rmtree(OUTDIR, ignore_errors=True)
for directory in (FIGDIR, FIGPRESENTDIR, FIGPAPERDIR, FIGSUPPDIR, PAPERTABLEDIR, SUPPTABLEDIR):
    os.makedirs(directory, exist_ok=True)

def sstr(x):
    return "" if pd.isna(x) else str(x)

def save_table_to(df, directory, name, index=False):
    os.makedirs(directory, exist_ok=True)
    path = os.path.join(directory, name)
    df.to_csv(path, index=index)
    return path

def tidy_round(df, digits=3):
    out = df.copy()
    for c in out.columns:
        if pd.api.types.is_numeric_dtype(out[c]):
            out[c] = out[c].round(digits)
    return out

def _excel_col(n):
    s = ""
    while n:
        n, r = divmod(n - 1, 26)
        s = chr(65 + r) + s
    return s

def _write_matrix_formula_table(ws, start_row, start_col, row_labels, col_labels, formula_builder, title, title_fmt, header_fmt, num_fmt):
    ws.write(start_row, start_col, title, title_fmt)
    ws.write(start_row + 1, start_col, "", header_fmt)
    for j, col in enumerate(col_labels):
        ws.write(start_row + 1, start_col + 1 + j, col, header_fmt)
    for i, row in enumerate(row_labels):
        ws.write(start_row + 2 + i, start_col, row, header_fmt)
        for j, col in enumerate(col_labels):
            ws.write_formula(start_row + 2 + i, start_col + 1 + j, formula_builder(row, col), num_fmt)

def export_excel_playground(df_source, output_path):
    if not EXPORT_EXCEL_PLAYGROUND:
        return None

    try:
        __import__("xlsxwriter")
    except Exception:
        print("Excel playground skipped: install xlsxwriter to create the .xlsx workbook.")
        return None

    excel_df = df_source.loc[df_source["possibility_flag"].eq(1), [
        "interaction_id", "type_simple", "E2_eff", "substituent", "conformer"
    ]].copy()
    excel_df = excel_df.rename(columns={"type_simple": "interaction_type", "E2_eff": "E2"})
    excel_df["E2"] = pd.to_numeric(excel_df["E2"], errors="coerce").fillna(0)
    excel_df = excel_df.sort_values(["substituent", "conformer", "interaction_type", "interaction_id"])

    conformers = [c for c in CONFORMER_ORDER if c in set(excel_df["conformer"])]
    substituents = [s for s in SUB_ORDER if s in set(excel_df["substituent"])]
    interaction_types = type_order(excel_df["interaction_type"])
    interaction_ids = sorted(excel_df["interaction_id"].dropna().astype(str).unique())

    with pd.ExcelWriter(output_path, engine="xlsxwriter") as writer:
        excel_df.to_excel(writer, sheet_name="Data", index=False)
        workbook = writer.book
        ws_data = writer.sheets["Data"]
        ws_summary = workbook.add_worksheet("Summary")
        ws_lists = workbook.add_worksheet("Lists")

        title_fmt = workbook.add_format({"bold": True, "font_size": 13, "font_color": "#FFFFFF", "bg_color": "#1F4E78", "align": "center"})
        header_fmt = workbook.add_format({"bold": True, "font_color": "#FFFFFF", "bg_color": "#5B9BD5", "border": 1, "align": "center"})
        subheader_fmt = workbook.add_format({"bold": True, "bg_color": "#D9EAF7", "border": 1})
        num_fmt = workbook.add_format({"num_format": "0.000", "border": 1})
        text_fmt = workbook.add_format({"border": 1})

        nrows = len(excel_df) + 1
        ncols = len(excel_df.columns)
        ws_data.add_table(0, 0, nrows - 1, ncols - 1, {
            "name": "EnergyData",
            "columns": [{"header": c} for c in excel_df.columns],
            "style": "Table Style Medium 2",
        })
        ws_data.freeze_panes(1, 0)
        ws_data.set_column("A:A", 42)
        ws_data.set_column("B:B", 16)
        ws_data.set_column("C:C", 12)
        ws_data.set_column("D:E", 14)

        list_blocks = [
            ("Conformers", conformers),
            ("InteractionTypes", interaction_types),
            ("Substituents", substituents),
            ("InteractionIDs", interaction_ids),
        ]
        for col_idx, (name, values) in enumerate(list_blocks):
            ws_lists.write(0, col_idx, name, header_fmt)
            for row_idx, value in enumerate(values, start=1):
                ws_lists.write(row_idx, col_idx, value)
        ws_lists.hide()

        if conformers:
            workbook.define_name("ConformerList", f"=Lists!$A$2:$A${len(conformers)+1}")
        if interaction_types:
            workbook.define_name("InteractionTypeList", f"=Lists!$B$2:$B${len(interaction_types)+1}")
        if substituents:
            workbook.define_name("SubstituentList", f"=Lists!$C$2:$C${len(substituents)+1}")
        if interaction_ids:
            workbook.define_name("InteractionIDList", f"=Lists!$D$2:$D${len(interaction_ids)+1}")

        dv_common = {
        "validate": "list",
        "ignore_blank": True,
        "show_error": False,
        }
        ws_summary.merge_range("A1:H1", "NCI Energy - Tabular Summary", title_fmt)
        ws_summary.set_column("A:A", 22)
        ws_summary.set_column("B:H", 15)

        controls = [
            ("Selected conformer", conformers[0] if conformers else ""),
            ("Selected interaction type", interaction_types[0] if interaction_types else ""),
            ("Selected substituent", substituents[0] if substituents else ""),
            ("Selected interaction id", interaction_ids[0] if interaction_ids else ""),
        ]
        ws_summary.write("A4", "Controls", title_fmt)
        for i, (label, value) in enumerate(controls, start=5):
            ws_summary.write(i-1, 0, label, subheader_fmt)
            ws_summary.write(i-1, 1, value, text_fmt)
        ws_summary.data_validation("B5", {**dv_common, "source": "=ConformerList"})
        ws_summary.data_validation("B6", {**dv_common, "source": "=InteractionTypeList"})
        ws_summary.data_validation("B7", {**dv_common, "source": "=SubstituentList"})
        ws_summary.data_validation("B8", {**dv_common, "source": "=InteractionIDList"})

        def _multi_match(table_col, control_cell):
            return (
                f'--ISNUMBER(SEARCH(","&EnergyData[{table_col}]&",",'
                f'","&SUBSTITUTE(SUBSTITUTE({control_cell},", ",",")," ,",",")&","))'
            )
        m_conf = _multi_match("conformer", "$B$5")
        m_type = _multi_match("interaction_type", "$B$6")
        m_sub = _multi_match("substituent", "$B$7")
        m_id = _multi_match("interaction_id", "$B$8")

        kpi_labels = [
            "Sum E2 for selected conformer",
            "Average E2 for selected conformer",
            "Sum E2 for conformer + type",
            "Average E2 for conformer + type",
            "Sum E2 for conformer + type + substituent",
            "Selected interaction E2 in selected substituent/conformer",
        ]
        kpi_formulas = [
            f'=SUMPRODUCT(EnergyData[E2]*{m_conf})',
            f'=IFERROR(SUMPRODUCT(EnergyData[E2]*{m_conf})/SUMPRODUCT({m_conf}),"")',
            f'=SUMPRODUCT(EnergyData[E2]*{m_conf}*{m_type})',
            f'=IFERROR(SUMPRODUCT(EnergyData[E2]*{m_conf}*{m_type})/SUMPRODUCT({m_conf}*{m_type}),"")',
            f'=SUMPRODUCT(EnergyData[E2]*{m_conf}*{m_type}*{m_sub})',
            f'=SUMPRODUCT(EnergyData[E2]*{m_id}*{m_sub}*{m_conf})',
        ]
        ws_summary.write("D4", "Quick formulas", title_fmt)
        for idx_k, (label, formula) in enumerate(zip(kpi_labels, kpi_formulas), start=5):
            ws_summary.write(idx_k-1, 3, label, subheader_fmt)
            ws_summary.write_formula(idx_k-1, 4, formula, num_fmt)

        row0 = 11
        _write_matrix_formula_table(
            ws_summary, row0, 0, substituents, conformers,
            lambda s, c: f'=IFERROR(AVERAGEIFS(EnergyData[E2],EnergyData[substituent],$A{row0+3+substituents.index(s)},EnergyData[conformer],{_excel_col(2+conformers.index(c))}${row0+2}),"")',
            "Average E2 by substituent × conformer", title_fmt, header_fmt, num_fmt
        )

        row1 = row0 + len(substituents) + 5
        _write_matrix_formula_table(
            ws_summary, row1, 0, substituents, conformers,
            lambda s, c: f'=SUMIFS(EnergyData[E2],EnergyData[substituent],$A{row1+3+substituents.index(s)},EnergyData[conformer],{_excel_col(2+conformers.index(c))}${row1+2})',
            "Sum E2 by substituent × conformer", title_fmt, header_fmt, num_fmt
        )

        row2 = row1 + len(substituents) + 5
        _write_matrix_formula_table(
            ws_summary, row2, 0, interaction_types, substituents,
            lambda t, s: f'=SUMIFS(EnergyData[E2],EnergyData[interaction_type],$A{row2+3+interaction_types.index(t)},EnergyData[substituent],{_excel_col(2+substituents.index(s))}${row2+2})',
            "Sum E2 by interaction type × substituent", title_fmt, header_fmt, num_fmt
        )

        row3 = row2 + len(interaction_types) + 5
        _write_matrix_formula_table(
            ws_summary, row3, 0, interaction_types, conformers,
            lambda t, c: f'=IFERROR(AVERAGEIFS(EnergyData[E2],EnergyData[interaction_type],$A{row3+3+interaction_types.index(t)},EnergyData[conformer],{_excel_col(2+conformers.index(c))}${row3+2}),"")',
            "Average E2 by interaction type × conformer", title_fmt, header_fmt, num_fmt
        )

        row4 = row3 + len(interaction_types) + 5
        _write_matrix_formula_table(
            ws_summary, row4, 0, substituents, conformers,
            lambda s, c: (
                f'=SUMPRODUCT(EnergyData[E2]*{m_id}'
                f'*--(EnergyData[substituent]=$A{row4+3+substituents.index(s)})'
                f'*--(EnergyData[conformer]={_excel_col(2+conformers.index(c))}${row4+2}))'
            ),
            "Selected interaction id: E2 by substituent × conformer", title_fmt, header_fmt, num_fmt
        )

        chart = workbook.add_chart({"type": "line"})
        for j, c in enumerate(conformers):
            chart.add_series({
                "name": ["Summary", row4 + 1, 1 + j],
                "categories": ["Summary", row4 + 2, 0, row4 + 1 + len(substituents), 0],
                "values": ["Summary", row4 + 2, 1 + j, row4 + 1 + len(substituents), 1 + j],
            })
        chart.set_title({"name": "Selected interaction profile"})
        chart.set_x_axis({"name": "Substituent"})
        chart.set_y_axis({"name": "E2"})
        chart.set_legend({"position": "bottom"})
        ws_summary.insert_chart(row4, max(7, len(conformers) + 3), chart, {"x_scale": 1.2, "y_scale": 1.1})

    print(f"Excel playground written: {output_path}")
    return output_path

def savefig(name):
    path = os.path.join(FIGDIR, name)
    plt.savefig(path, dpi=200)
    plt.close()
    return path

def savefig_presentation(name):
    path = os.path.join(FIGPRESENTDIR, name)
    plt.savefig(path, dpi=230)
    plt.close()
    return path

def savefig_paper(name):
    path = os.path.join(FIGPAPERDIR, name)
    plt.savefig(path, dpi=1200, bbox_inches="tight", facecolor=plt.gcf().get_facecolor())
    plt.close()
    return path

def savefig_supplementary(name):
    path = os.path.join(FIGSUPPDIR, name)
    plt.savefig(path, dpi=1200, bbox_inches="tight", facecolor=plt.gcf().get_facecolor())
    plt.close()
    return path

def norm_type(x):
    s = sstr(x).lower()
    s = s.replace("π", "pi").replace("σ", "sigma").replace("∗", "*")
    s = s.replace("→", "-").replace("->", "-").replace(" ", "")
    return s

def _orbital_family_token(token):
    s = norm_type(token)
    s = s.replace("*", "").replace("star", "")
    if s.startswith("n") or "lp" in s:
        return "n"
    if "pi" in s:
        return "π"
    if "sigma" in s or "sig" in s:
        return "σ"
    return None

def type_simple(x):
    raw = sstr(x)
    s = norm_type(raw)

    parts = [p for p in re.split(r"-+", s) if p != ""]
    if len(parts) >= 2:
        donor = _orbital_family_token(parts[0])
        acceptor = _orbital_family_token(parts[-1])
        has_acceptor_star = "*" in parts[-1] or "star" in parts[-1]
        if donor and acceptor and has_acceptor_star:
            return f"{donor}→{acceptor}*"

    donor_n = s.startswith("n") or s.startswith("lp")
    has_pi = "pi" in s
    has_sig = "sigma" in s or "sig" in s
    has_star = "*" in s or "star" in s
    if donor_n and has_pi and has_star:
        return "n→π*"
    if donor_n and has_sig and has_star:
        return "n→σ*"
    if has_sig and has_pi and has_star:
        return "σ→π*"
    if has_sig and has_star:
        return "σ→σ*"
    return raw

def conformer_family(conf):
    c = sstr(conf)
    if c.startswith("P"):
        return "P"
    if c.startswith("N"):
        return "N"
    return "other"

def ensure_interaction_id(df):
    if "interaction_id" in df.columns:
        return df
    cols = ["donor_orbital", "acceptor_orbital", "donor_labels", "acceptor_labels", "interaction_type"]
    if all(c in df.columns for c in cols):
        df["interaction_id"] = (
            df["donor_orbital"].astype(str) + "__" +
            df["acceptor_orbital"].astype(str) + "__" +
            df["donor_labels"].astype(str) + "__" +
            df["acceptor_labels"].astype(str) + "__" +
            df["interaction_type"].astype(str)
        )
    else:
        df["interaction_id"] = (
            df["donor_labels"].astype(str) + "__" +
            df["acceptor_labels"].astype(str) + "__" +
            df["interaction_type"].astype(str)
        )
    return df

def long_label(row):
    if "donor_orbital" in row.index and "acceptor_orbital" in row.index:
        return f"{row.get('donor_orbital','')} → {row.get('acceptor_orbital','')} | {row.get('donor_labels','')} → {row.get('acceptor_labels','')} | {row.get('interaction_type','')}"
    return f"{row.get('donor_labels','')} → {row.get('acceptor_labels','')} | {row.get('interaction_type','')}"

def corr_pair(x, y, min_n=MIN_CORR_N):
    x = pd.Series(x)
    y = pd.Series(y)
    mask = x.notna() & y.notna()
    n = int(mask.sum())
    if n < min_n:
        return np.nan, n
    a = x[mask].astype(float)
    b = y[mask].astype(float)
    if a.std(ddof=0) == 0 or b.std(ddof=0) == 0:
        return np.nan, n
    return float(np.corrcoef(a, b)[0, 1]), n

def type_order(vals):
    pref = ["n→π*", "n→σ*", "π→σ*", "σ→π*", "σ→σ*"]
    present = list(pd.Series(vals).dropna().unique())
    return [x for x in pref if x in present] + [x for x in sorted(present) if x not in pref]

def atom_numbers_from_label(label):
    return [int(n) for n in re.findall(r"[A-Z][a-z]?(\d+)", sstr(label))]

def is_core_interaction(row):
    nums = atom_numbers_from_label(row.get("donor_labels", "")) + atom_numbers_from_label(row.get("acceptor_labels", ""))
    return bool(nums) and all(n <= CORE_MAX_ATOM for n in nums)

def is_internal_1_9_only_label(label):
    nums = atom_numbers_from_label(label)
    return bool(nums) and all(1 <= n <= 9 for n in nums)

def has_internal_1_9_only_label(row):
    return (
        is_internal_1_9_only_label(row.get("donor_labels", ""))
        or is_internal_1_9_only_label(row.get("acceptor_labels", ""))
    )

def add_hammett_top_axis(ax, ordered_substituents, show_label=True):
    try:
        secax = ax.secondary_xaxis("top")
        ticks = np.arange(len(ordered_substituents))
        secax.set_xticks(ticks)
        secax.set_xticklabels([f"{HAMMETT_SIGMA_P.get(s, np.nan):.2f}" for s in ordered_substituents], rotation=45, ha="left", fontsize=8)
        if show_label:
            secax.set_xlabel("Hammett σp")
    except Exception:
        pass

STRUCTURE_CACHE = {}
BOND_LENGTH_SAME_TOL_A = 0.03
COVALENT_RADII = {
    "H": 0.37, "C": 0.77, "N": 0.75, "O": 0.73,
    "F": 0.71, "P": 1.07, "S": 1.05, "Cl": 1.02,
    "Br": 1.20, "I": 1.39, "Si": 1.11, "B": 0.85,
    "Na": 1.66, "K": 2.03, "Ca": 1.74, "Fe": 1.26,
    "Zn": 1.22, "Cu": 1.32, "Mg": 1.41, "Al": 1.18,
}
MAX_VALENCE = {"H": 1, "F": 1, "Cl": 1, "Br": 1, "I": 1, "O": 2, "N": 4, "C": 4, "B": 3, "P": 5, "S": 6}
LITERATURE_BOND_LENGTHS = {
    ("C", "H", "single"): 1.09, ("N", "H", "single"): 1.01, ("O", "H", "single"): 0.96,
    ("C", "C", "single"): 1.54, ("C", "C", "double"): 1.34,
    ("C", "N", "single"): 1.47, ("C", "N", "double"): 1.28,
    ("C", "O", "single"): 1.43, ("C", "O", "double"): 1.21,
    ("C", "F", "single"): 1.35, ("C", "Cl", "single"): 1.77,
}

def parse_xyz_structure(path):
    with open(path, "r", encoding="utf-8") as fh:
        lines = [ln.split() for ln in fh if ln.strip()]
    start = 2 if lines and lines[0][0].isdigit() else 0
    return [(p[0], float(p[1]), float(p[2]), float(p[3])) for p in lines[start:] if len(p) >= 4]

def parse_pdb_structure(path):
    records, serial_to_idx, serial_bonds = [], {}, {}
    with open(path, "r", encoding="utf-8") as fh:
        for ln in fh:
            if ln.startswith(("ATOM", "HETATM")):
                try:
                    serial = int(ln[6:11])
                    name = ln[12:16].strip()
                    elem = (ln[76:78].strip() if len(ln) >= 78 else "") or re.sub(r"[^A-Za-z]", "", name)[:2].title()
                    if elem not in COVALENT_RADII and elem[:1] in COVALENT_RADII:
                        elem = elem[:1]
                    records.append((elem, float(ln[30:38]), float(ln[38:46]), float(ln[46:54])))
                    serial_to_idx[serial] = len(records)
                except Exception:
                    continue
            elif ln.startswith("CONECT"):
                nums = [int(n) for n in ln.split()[1:] if n.isdigit()]
                if nums:
                    src = nums[0]
                    for tgt in nums[1:]:
                        serial_bonds.setdefault(src, set()).add(tgt)
                        serial_bonds.setdefault(tgt, set()).add(src)
    adj = {i: set() for i in range(1, len(records) + 1)}
    for src, tgts in serial_bonds.items():
        if src in serial_to_idx:
            adj[serial_to_idx[src]].update(serial_to_idx[t] for t in tgts if t in serial_to_idx)
    return records, adj if any(adj[i] for i in adj) else None

def atom_label_for_index(records, idx):
    return f"{records[idx - 1][0]}{idx}"

def atom_tokens(label):
    return re.findall(r"[A-Z][a-z]?\d+", sstr(label))

def atom_index(atom_label):
    m = re.match(r"[A-Z][a-z]?(\d+)", sstr(atom_label))
    return int(m.group(1)) if m else None

def atom_element(atom_label):
    m = re.match(r"([A-Z][a-z]?)(\d+)", sstr(atom_label))
    return m.group(1) if m else None

def h_substituent_geometry_token(raw_token, substituent, records=None, structure_path=None):
    elem, idx = atom_element(raw_token), atom_index(raw_token)
    if idx is None or elem is None:
        return raw_token, ""

    is_h_xyz = sstr(substituent) == "H" and sstr(structure_path).lower().endswith(".xyz")
    if is_h_xyz and elem == "H" and idx == 28:
        geom_idx = 8
        rule = "special H28->xyz8"
    elif is_h_xyz and idx >= 8:
        geom_idx = idx + 1
        rule = "reverse +1 after corrected index 8"
    else:
        geom_idx = idx
        rule = ""

    geom_token = f"{elem}{geom_idx}"
    note = "" if geom_idx == idx else f"H mainframe-to-XYZ reverse mapping ({raw_token}->{geom_token}; {rule})"
    if records and not (1 <= geom_idx <= len(records)):
        note = f"{note}; mapped index outside structure range" if note else "index outside structure range"
    return geom_token, note

def display_token_from_geometry(geom_token, substituent, structure_path=None):
    elem, idx = atom_element(geom_token), atom_index(geom_token)
    is_h_xyz = sstr(substituent) == "H" and sstr(structure_path).lower().endswith(".xyz")
    if is_h_xyz and elem and idx:
        if elem == "H" and idx == 8:
            return "H28"
        if idx > 8:
            return f"{elem}{idx - 1}"
    return geom_token

def atom_distance(records, i, j):
    _, xi, yi, zi = records[i - 1]
    _, xj, yj, zj = records[j - 1]
    return float(math.dist((xi, yi, zi), (xj, yj, zj)))

def build_bond_graph(records, explicit_adj=None, add_tol=0.42):
    n = len(records)
    if explicit_adj and any(explicit_adj.get(i) for i in explicit_adj):
        return {i: set(explicit_adj.get(i, set())) for i in range(1, n + 1)}, "PDB_CONECT"
    adj, lengths = {i: set() for i in range(1, n + 1)}, {}
    for i, (el_i, xi, yi, zi) in enumerate(records, start=1):
        for j, (el_j, xj, yj, zj) in enumerate(records[i:], start=i + 1):
            if el_i == "H" and el_j == "H":
                continue
            dist = math.dist((xi, yi, zi), (xj, yj, zj))
            if dist <= COVALENT_RADII.get(el_i, 0.77) + COVALENT_RADII.get(el_j, 0.77) + add_tol:
                adj[i].add(j); adj[j].add(i); lengths[tuple(sorted((i, j)))] = dist
    changed = True
    while changed:
        changed = False
        for i, (el_i, *_xyz) in enumerate(records, start=1):
            while len(adj[i]) > MAX_VALENCE.get(el_i, 4):
                j = max(adj[i], key=lambda k: lengths.get(tuple(sorted((i, k))), atom_distance(records, i, k)))
                adj[i].remove(j); adj[j].discard(i); changed = True
    return adj, "covalent_radii_pruned"

def find_structure_path(substituent, conformer):
    for folder in ["", "xyz", "structures", "coordinates"]:
        for base in [f"{substituent}_{conformer}", f"{conformer}_{substituent}"]:
            for ext in [".pdb", ".xyz"]:
                path = os.path.join(folder, base + ext) if folder else base + ext
                if os.path.exists(path):
                    return path
    return None

def read_structure_cached(substituent, conformer):
    key = (substituent, conformer)
    if key in STRUCTURE_CACHE:
        return STRUCTURE_CACHE[key]
    path = find_structure_path(substituent, conformer)
    if not path:
        out = {"path": None, "records": [], "coords": {}, "adj": {}, "graph_source": "missing_structure"}
    else:
        try:
            records, explicit_adj = parse_pdb_structure(path) if path.lower().endswith(".pdb") else (parse_xyz_structure(path), None)
            adj, graph_source = build_bond_graph(records, explicit_adj=explicit_adj)
            coords = {atom_label_for_index(records, i): np.array(records[i - 1][1:4], dtype=float) for i in range(1, len(records) + 1)}
            out = {"path": path, "records": records, "coords": coords, "adj": adj, "graph_source": graph_source}
        except Exception as exc:
            out = {"path": path, "records": [], "coords": {}, "adj": {}, "graph_source": f"parse_failed: {exc}"}
    STRUCTURE_CACHE[key] = out
    return out

def choose_neighbor_for_single_atom(records, adj, idx):
    if idx not in adj or not adj[idx]:
        return None, "no bonded neighbor in inferred graph"
    elem, neighbors = records[idx - 1][0], list(adj[idx])
    def score(j):
        ej = records[j - 1][0]
        return (
            0 if ej != "H" else 5,
            0 if not (1 <= idx <= 9) else (0 if 1 <= j <= 9 else 2),
            -1 if elem == "C" and ej == "C" else 0 if ej != "H" else 2,
            atom_distance(records, idx, j),
            j,
        )
    return min(neighbors, key=score), "single-label inferred from covalent graph: heavy/inward/shortest heuristic"

def resolve_bond_atoms_from_label(struct, label, row):
    records, adj, sub = struct["records"], struct["adj"], row.get("substituent", "")
    raw_tokens = atom_tokens(label)
    if not raw_tokens:
        return [], [], "no atom token in label", ""
    geom_tokens, notes = [], []
    for tok in raw_tokens[:2]:
        gt, note = h_substituent_geometry_token(tok, sub, records, struct.get("path"))
        geom_tokens.append(gt)
        if note:
            notes.append(note)
    note_txt = "; ".join(notes)
    if len(raw_tokens) >= 2:
        method = "explicit two-atom label" if len(raw_tokens) == 2 else "multi-token label; first two atoms used"
        return raw_tokens[:2], geom_tokens[:2], method, note_txt
    idx = atom_index(geom_tokens[0])
    if not idx or idx < 1 or idx > len(records):
        return [], [], "single atom not found in structure", note_txt
    partner, method = choose_neighbor_for_single_atom(records, adj, idx)
    if not partner:
        return [], [], method, note_txt
    partner_geom = atom_label_for_index(records, partner)
    return [raw_tokens[0], display_token_from_geometry(partner_geom, sub, struct.get("path"))], [geom_tokens[0], partner_geom], method, note_txt

def explicit_bond_length(struct, geom_atoms):
    records = struct["records"]
    idxs = [atom_index(a) for a in geom_atoms]
    if len(idxs) != 2 or any(i is None or i < 1 or i > len(records) for i in idxs):
        return np.nan
    return atom_distance(records, idxs[0], idxs[1])

def bond_order_hint(row, role, atoms):
    typ = sstr(row.get("type_simple", ""))
    fam = typ.split("→")[0] if role == "donor" and "→" in typ else typ.split("→")[-1].replace("*", "") if "→" in typ else ""
    elems = sorted(filter(None, [atom_element(a) for a in atoms]))
    pair = tuple(elems)
    if fam == "π":
        return "double"
    if pair == ("C", "O") and role == "acceptor" and typ in {"n→π*", "σ→π*"}:
        return "double"
    if pair == ("C", "C") and role == "acceptor" and typ == "σ→π*":
        return "double"
    if pair == ("C", "N") and role == "acceptor" and typ == "σ→π*":
        return "double"
    return "single"

def literature_bond_length(row, role, display_atoms):
    if len(display_atoms) != 2:
        return np.nan, "unresolved bond"
    elems = sorted(filter(None, [atom_element(a) for a in display_atoms]))
    if len(elems) != 2:
        return np.nan, "unparsed elements"
    order = bond_order_hint(row, role, display_atoms)
    ref = LITERATURE_BOND_LENGTHS.get((elems[0], elems[1], order)) or LITERATURE_BOND_LENGTHS.get((elems[0], elems[1], "single"))
    if ref is None:
        return np.nan, f"no built-in reference for {elems[0]}-{elems[1]}"
    used_order = order if (elems[0], elems[1], order) in LITERATURE_BOND_LENGTHS else "single"
    return float(ref), f"typical literature {elems[0]}-{elems[1]} {used_order} bond"

def classify_bond_length_change(delta):
    if pd.isna(delta):
        return "unresolved"
    if delta < -BOND_LENGTH_SAME_TOL_A:
        return "shorter"
    if delta > BOND_LENGTH_SAME_TOL_A:
        return "longer"
    return "same"

def bond_length_comment(status, delta, method):
    if status == "shorter":
        return f"shorter than literature by {abs(delta):.3f} Å; {method}"
    if status == "longer":
        return f"longer than literature by {abs(delta):.3f} Å; {method}"
    if status == "same":
        return f"approximately same as literature (within ±{BOND_LENGTH_SAME_TOL_A:.2f} Å); {method}"
    return f"not compared; {method}"

def build_bond_length_literature_table(df_in):
    work = df_in.copy()
    meta_cols = [c for c in ["interaction_id", "interaction_code", "donor_labels", "acceptor_labels", "interaction_type", "type_simple", "core_interaction", "has_internal_1_9_only_label"] if c in work.columns]
    interactions = work.sort_values(["interaction_id", "substituent", "conformer"]).drop_duplicates("interaction_id")[meta_cols].copy()
    sample_grid = pd.MultiIndex.from_product([SUB_ORDER, CONFORMER_ORDER], names=["substituent", "conformer"]).to_frame(index=False)
    full = sample_grid.assign(_key=1).merge(interactions.assign(_key=1), on="_key").drop(columns="_key")
    vals = work[["substituent", "conformer", "interaction_id", "exists_flag", "possibility_flag", "E2_kcal/mol"]].copy()
    vals["E2_kcal/mol"] = pd.to_numeric(vals["E2_kcal/mol"], errors="coerce").fillna(0)
    vals = vals.groupby(["substituent", "conformer", "interaction_id"], observed=False).agg(
        exists_flag=("exists_flag", "max"), possibility_flag=("possibility_flag", "max"), E2_kcal_mol=("E2_kcal/mol", "max")
    ).reset_index()
    full = full.merge(vals, on=["substituent", "conformer", "interaction_id"], how="left")
    full[["exists_flag", "possibility_flag"]] = full[["exists_flag", "possibility_flag"]].fillna(0).astype(int)
    full["E2_kcal_mol"] = full["E2_kcal_mol"].fillna(0)
    full["E2_for_table"] = np.where(full["exists_flag"].eq(1) & full["possibility_flag"].eq(1), full["E2_kcal_mol"], 0.0)

    recs = []
    for row_s in full.to_dict("records"):
        struct = read_structure_cached(row_s["substituent"], row_s["conformer"])
        for role, label_col in [("donor", "donor_labels"), ("acceptor", "acceptor_labels")]:
            label = row_s.get(label_col, "")
            disp_atoms, geom_atoms, method, correction_note = resolve_bond_atoms_from_label(struct, label, row_s)
            obs = explicit_bond_length(struct, geom_atoms) if len(geom_atoms) == 2 else np.nan
            ref, ref_note = literature_bond_length(row_s, role, disp_atoms)
            delta = obs - ref if pd.notna(obs) and pd.notna(ref) else np.nan
            status = classify_bond_length_change(delta)
            recs.append({
                "substituent": row_s.get("substituent", ""), "conformer": row_s.get("conformer", ""),
                "interaction_code": row_s.get("interaction_code", ""), "interaction_id": row_s.get("interaction_id", ""),
                "role": role.title(), "original_label": label, "resolved_bond_label": "|".join(disp_atoms) if len(disp_atoms) == 2 else "",
                "geometry_bond_label_used_for_distance": "|".join(geom_atoms) if len(geom_atoms) == 2 else "",
                "elements": "-".join(sorted(filter(None, [atom_element(a) for a in disp_atoms]))) if len(disp_atoms) == 2 else "",
                "interaction_type": row_s.get("type_simple", ""), "exists_flag": row_s.get("exists_flag", 0),
                "possibility_flag": row_s.get("possibility_flag", 0), "E2_for_table": row_s.get("E2_for_table", 0.0),
                "core_interaction": row_s.get("core_interaction", False), "has_internal_1_9_only_label": row_s.get("has_internal_1_9_only_label", False),
                "structure_file": struct.get("path"), "bond_graph_source": struct.get("graph_source"), "resolution_method": method,
                "atom_label_correction_note": correction_note, "observed_bond_length_A": obs, "literature_bond_length_A": ref,
                "delta_obs_minus_lit_A": delta, "bond_length_status": status,
                "bond_length_comment": bond_length_comment(status, delta, method), "reference_note": ref_note,
            })
    out = pd.DataFrame(recs)
    full_expected = len(out)
    out = out[out["exists_flag"].eq(1) | out["bond_length_status"].eq("same")].copy()
    out.attrs["full_grid_expected_entries"] = full_expected
    return out[[
        "substituent", "conformer", "interaction_code", "interaction_id", "role", "original_label", "resolved_bond_label",
        "geometry_bond_label_used_for_distance", "elements",
        "interaction_type", "exists_flag", "possibility_flag", "E2_for_table", "core_interaction", "has_internal_1_9_only_label",
        "structure_file", "bond_graph_source", "resolution_method", "atom_label_correction_note", "observed_bond_length_A", "literature_bond_length_A",
        "delta_obs_minus_lit_A", "bond_length_status", "bond_length_comment", "reference_note",
    ]]

def _bond_summary_group(bond_table, group_cols):
    if bond_table is None or bond_table.empty:
        return pd.DataFrame({"note": ["No bond-length rows were generated."]})

    value_cols = dict(
        n_entries=("bond_length_status", "size"),
        n_existing=("exists_flag", "sum"),
        mean_delta_A=("delta_obs_minus_lit_A", "mean"),
        mean_observed_A=("observed_bond_length_A", "mean"),
        mean_literature_A=("literature_bond_length_A", "mean"),
    )

    if group_cols:
        base = bond_table.groupby(group_cols, observed=False).agg(**value_cols).reset_index()
        counts = (
            bond_table.pivot_table(
                index=group_cols,
                columns="bond_length_status",
                values="interaction_id",
                aggfunc="count",
                fill_value=0,
                observed=False,
            )
            .reset_index()
        )
        out = base.merge(counts, on=group_cols, how="left")
    else:
        out = pd.DataFrame([{
            "n_entries": len(bond_table),
            "n_existing": bond_table["exists_flag"].sum(),
            "mean_delta_A": bond_table["delta_obs_minus_lit_A"].mean(),
            "mean_observed_A": bond_table["observed_bond_length_A"].mean(),
            "mean_literature_A": bond_table["literature_bond_length_A"].mean(),
        }])
        counts = bond_table["bond_length_status"].value_counts().to_dict()
        for status, n in counts.items():
            out[status] = int(n)

    for col in ["shorter", "longer", "same", "unresolved"]:
        if col not in out.columns:
            out[col] = 0
        out[f"fraction_{col}"] = out[col] / out["n_entries"].replace(0, np.nan)
    return out

def build_bond_length_literature_summary(bond_table):
    if bond_table is None or bond_table.empty:
        return pd.DataFrame({"note": ["No bond-length rows were generated."]})

    meta = pd.DataFrame([{
        "summary_level": "row_count_check",
        "bond_length_status": "filtered_expected_vs_observed",
        "n_entries": len(bond_table),
        "full_grid_expected_entries_before_filter": bond_table.attrs.get("full_grid_expected_entries", np.nan),
        "filter_rule": "kept exists_flag=1; additionally kept exists_flag=0 only when bond_length_status == same",
        "n_unique_interactions_included": bond_table["interaction_id"].nunique(),
        "n_substituents_included": bond_table["substituent"].nunique(),
        "n_conformers_included": bond_table["conformer"].nunique(),
    }])

    overall = _bond_summary_group(bond_table, [])
    overall.insert(0, "summary_level", "overall")

    type_role = _bond_summary_group(bond_table, ["interaction_type", "role"])
    type_role.insert(0, "summary_level", "interaction_type_role")

    return pd.concat([meta, overall, type_role], ignore_index=True, sort=False)

def build_bond_length_extra_summaries(bond_table):
    return {
        "Table_bond_length_summary_by_substituent_conformer.csv": _bond_summary_group(bond_table, ["substituent", "conformer"]),
        "Table_bond_length_top_deviations.csv": bond_table.dropna(subset=["delta_obs_minus_lit_A"]).assign(abs_delta_A=lambda x: x["delta_obs_minus_lit_A"].abs()).sort_values("abs_delta_A", ascending=False).head(80),
    }

def build_summary_table(df_in):
    summary_local = df_in.groupby(["substituent", "conformer"], observed=False).agg(
        sigma_p=("sigma_p", "first"),
        conformer_family=("conformer_family", "first"),
        possible_interactions=("possibility_flag", "sum"),
        existing_interactions=("exists_flag", "sum"),
        total_E2=("E2_eff", "sum"),
        mean_E2_existing=("E2_present", "mean"),
        max_E2=("E2_present", "max"),
    ).reset_index()

    npi_tmp_local = df_in[df_in["is_npi"]].groupby(["substituent", "conformer"], observed=False).agg(
        npi_total_E2=("E2_eff", "sum"),
        npi_count=("exists_flag", "sum"),
        npi_max_E2=("E2_present", "max"),
        npi_mean_distance=("distance", "mean") if "distance" in df_in.columns else ("E2_eff", "mean"),
        npi_mean_angle=("angle", "mean") if "angle" in df_in.columns else ("E2_eff", "mean"),
        npi_mean_dihedral=("dihedral", "mean") if "dihedral" in df_in.columns else ("E2_eff", "mean"),
    ).reset_index()

    summary_local = summary_local.merge(npi_tmp_local, on=["substituent", "conformer"], how="left")
    for c in ["npi_total_E2", "npi_count", "npi_max_E2"]:
        summary_local[c] = summary_local[c].fillna(0)
    summary_local["non_npi_total_E2"] = summary_local["total_E2"] - summary_local["npi_total_E2"]
    summary_local["npi_fraction_E2"] = summary_local["npi_total_E2"] / summary_local["total_E2"].replace(0, np.nan)
    summary_local["existence_fraction"] = summary_local["existing_interactions"] / summary_local["possible_interactions"].replace(0, np.nan)
    return summary_local.sort_values(["sigma_p", "conformer"])

# fixed csv layout, no need to be clever here
df = pd.read_csv(FILE)
df = ensure_interaction_id(df)

for c in ["E2_kcal/mol", "distance", "angle", "dihedral"]:
    if c in df.columns:
        df[c] = pd.to_numeric(df[c], errors="coerce")

df["exists_flag"] = pd.to_numeric(df["exists_flag"], errors="coerce").fillna(0).astype(int)
df["possibility_flag"] = pd.to_numeric(df["possibility_flag"], errors="coerce").fillna(0).astype(int)
df["substituent"] = df["substituent"].astype(str)
df["conformer"] = df["conformer"].astype(str)
df["type_simple"] = df["interaction_type"].apply(type_simple)
EXPECTED_MAJOR_TYPES = ["n→π*", "n→σ*", "π→σ*", "σ→π*", "σ→σ*"]

df["is_npi"] = df["type_simple"].eq("n→π*")
df["core_interaction"] = df.apply(is_core_interaction, axis=1)
df["has_internal_1_9_only_label"] = df.apply(has_internal_1_9_only_label, axis=1)
df["conformer_family"] = df["conformer"].apply(conformer_family)
df["sigma_p"] = df["substituent"].map(HAMMETT_SIGMA_P)
df["E2_eff"] = np.where(df["possibility_flag"].eq(1), df["E2_kcal/mol"].fillna(0), np.nan)
df["E2_present"] = np.where(df["exists_flag"].eq(1) & df["possibility_flag"].eq(1), df["E2_kcal/mol"], np.nan)
df["sub_order"] = df["substituent"].map({s: i for i, s in enumerate(SUB_ORDER)})
df["conf_order"] = df["conformer"].map({c: i for i, c in enumerate(CONFORMER_ORDER)})
df = df.sort_values(["sub_order", "conf_order", "type_simple", "interaction_id"])

def _prepare_secondary_bubble_df(file_path):
    x = pd.read_csv(file_path)
    x = ensure_interaction_id(x)
    for c in ["E2_kcal/mol", "distance", "angle", "dihedral"]:
        x[c] = pd.to_numeric(x[c], errors="coerce")
    x["exists_flag"] = 1
    x["possibility_flag"] = 1
    x["conformer"] = "Biomotifs"
    x["interaction_type"] = x["interaction_type"].fillna("").astype(str)
    x["donor_labels"] = x["donor_labels"].fillna("").astype(str)
    x["acceptor_labels"] = x["acceptor_labels"].fillna("").astype(str)
    x["type_simple"] = x["interaction_type"].apply(type_simple)
    x["conformer_family"] = "other"
    x["E2_eff"] = x["E2_kcal/mol"].fillna(0)
    x["E2_present"] = x["E2_kcal/mol"]
    x["core_interaction"] = False
    return x

# motif rows only feed the geometry plots below
df_biomotifs = _prepare_secondary_bubble_df(MOTIF_FILE)

base = df.drop_duplicates("interaction_id").copy()
base["long_label"] = base.apply(long_label, axis=1)
contrib = df.groupby("interaction_id").agg(total_E2=("E2_eff", "sum")).reset_index()
dictionary = base[["interaction_id", "long_label", "type_simple", "is_npi", "core_interaction"]].merge(contrib, on="interaction_id", how="left")
dictionary = dictionary.sort_values(["is_npi", "type_simple", "total_E2"], ascending=[False, True, False]).reset_index(drop=True)
dictionary["interaction_code"] = [f"I{i+1:02d}" for i in range(len(dictionary))]
code_map = dictionary.set_index("interaction_id")["interaction_code"].to_dict()
label_map = dictionary.set_index("interaction_id")["long_label"].to_dict()
type_map = dictionary.set_index("interaction_id")["type_simple"].to_dict()
df["interaction_code"] = df["interaction_id"].map(code_map)

summary = build_summary_table(df)

summary_core = build_summary_table(df[df["core_interaction"]].copy())

bond_lit_table = build_bond_length_literature_table(df)
bond_lit_path = os.path.join(PAPERTABLEDIR, "Table_bond_length_vs_literature.csv")
tidy_round(bond_lit_table).to_csv(bond_lit_path, index=False)

bond_lit_summary = build_bond_length_literature_summary(bond_lit_table)
bond_lit_summary_path = os.path.join(PAPERTABLEDIR, "Table_bond_length_vs_literature_summary.csv")
tidy_round(bond_lit_summary).to_csv(bond_lit_summary_path, index=False)

for _bond_name, _bond_df in build_bond_length_extra_summaries(bond_lit_table).items():
    tidy_round(_bond_df).to_csv(os.path.join(PAPERTABLEDIR, _bond_name), index=False)

h_xyz_audit = bond_lit_table[
    bond_lit_table["substituent"].eq("H")
    & bond_lit_table["structure_file"].astype(str).str.lower().str.endswith(".xyz")
    & bond_lit_table["atom_label_correction_note"].astype(str).ne("")
][[
    "substituent", "conformer", "role", "original_label",
    "resolved_bond_label", "geometry_bond_label_used_for_distance",
    "atom_label_correction_note", "observed_bond_length_A",
    "literature_bond_length_A", "delta_obs_minus_lit_A",
    "bond_length_status"
]].copy()
tidy_round(h_xyz_audit).to_csv(os.path.join(PAPERTABLEDIR, "Table_H_xyz_reverse_label_audit.csv"), index=False)

type_sample = df[df["possibility_flag"].eq(1)].groupby(["substituent", "sigma_p", "conformer", "conformer_family", "type_simple"], observed=False).agg(total_E2=("E2_eff", "sum"), count_existing=("exists_flag", "sum"), mean_E2=("E2_present", "mean"), max_E2=("E2_present", "max")).reset_index()

type_sample_core = df[df["possibility_flag"].eq(1) & df["core_interaction"]].groupby(["substituent", "sigma_p", "conformer", "conformer_family", "type_simple"], observed=False).agg(total_E2=("E2_eff", "sum"), count_existing=("exists_flag", "sum"), mean_E2=("E2_present", "mean"), max_E2=("E2_present", "max")).reset_index()

def _safe_sheet_name(name, used):
    clean = re.sub(r"[\\/*?:\[\]]", "_", str(name))[:31] or "Sheet"
    base = clean
    i = 1
    while clean in used:
        suffix = f"_{i}"
        clean = base[:31-len(suffix)] + suffix
        i += 1
    used.add(clean)
    return clean

def _paper_formats(wb):
    return {
        "title": wb.add_format({"bold": True, "font_size": 14, "font_color": "#FFFFFF", "bg_color": "#1F4E78", "align": "center", "valign": "vcenter"}),
        "section": wb.add_format({"bold": True, "font_size": 11, "font_color": "#FFFFFF", "bg_color": "#305496", "align": "left", "valign": "vcenter"}),
        "note": wb.add_format({"italic": True, "font_color": "#666666", "text_wrap": True}),
        "header": wb.add_format({"bold": True, "font_color": "#FFFFFF", "bg_color": "#5B9BD5", "border": 1, "align": "center", "valign": "vcenter"}),
        "label": wb.add_format({"bold": True, "bg_color": "#D9EAF7", "border": 1, "align": "center", "valign": "vcenter"}),
        "num": wb.add_format({"num_format": "0.000", "border": 1, "align": "center"}),
        "text": wb.add_format({"border": 1, "text_wrap": True, "valign": "top"}),
        "blank": wb.add_format({"border": 1}),
    }

def _table_source(df_in):
    work = df_in.copy()
    meta_cols = [
        "interaction_id", "type_simple", "interaction_type", "donor_labels", "acceptor_labels",
        "core_interaction", "has_internal_1_9_only_label"
    ]
    meta_cols = [c for c in meta_cols if c in work.columns]
    interaction_meta = work.sort_values("interaction_id").drop_duplicates("interaction_id")[meta_cols].copy()

    sample_grid = pd.MultiIndex.from_product(
        [SUB_ORDER, CONFORMER_ORDER], names=["substituent", "conformer"]
    ).to_frame(index=False)
    sample_grid["sigma_p"] = sample_grid["substituent"].map(HAMMETT_SIGMA_P)
    sample_grid["conformer_family"] = sample_grid["conformer"].apply(conformer_family)

    grid = sample_grid.assign(_key=1).merge(interaction_meta.assign(_key=1), on="_key").drop(columns="_key")

    val = work[["substituent", "conformer", "interaction_id", "exists_flag", "possibility_flag", "E2_kcal/mol"]].copy()
    val["E2_kcal/mol"] = pd.to_numeric(val["E2_kcal/mol"], errors="coerce").fillna(0)
    val = val.groupby(["substituent", "conformer", "interaction_id"], observed=False).agg(
        exists_flag=("exists_flag", "max"),
        possibility_flag=("possibility_flag", "max"),
        E2_raw=("E2_kcal/mol", "max"),
    ).reset_index()

    out = grid.merge(val, on=["substituent", "conformer", "interaction_id"], how="left")
    out["exists_flag"] = out["exists_flag"].fillna(0).astype(int)
    out["possibility_flag"] = out["possibility_flag"].fillna(0).astype(int)
    out["E2_raw"] = out["E2_raw"].fillna(0)
    active = out["exists_flag"].eq(1) & out["possibility_flag"].eq(1)
    out["E2_table"] = np.where(active, out["E2_raw"], 0.0)
    out["count_table"] = np.where(active, 1, 0)
    return out

def _scope_df(core_only=False, partition=None):
    work = _table_source(df)
    if core_only:
        work = work[work["core_interaction"].eq(True)].copy()
    if partition == "ring":
        work = work[work["has_internal_1_9_only_label"].eq(True)].copy()
    elif partition == "nonring":
        work = work[work["has_internal_1_9_only_label"].eq(False)].copy()
    return work

def _sum_matrix(work, typ, index_col, index_order, index_name):
    g = work[work["type_simple"].eq(typ)].copy()
    table = g.pivot_table(index=index_col, columns="substituent", values="E2_table", aggfunc="sum", fill_value=0)
    table = table.reindex(index=index_order, columns=SUB_ORDER).fillna(0)
    return table.reset_index().rename(columns={index_col: index_name})

def _family_matrix(work, typ):
    return _sum_matrix(work, typ, "conformer_family", ["N", "P"], "Conformer Family")

def _conformer_matrix(work, typ):
    return _sum_matrix(work, typ, "conformer", CONFORMER_ORDER, "Conformer")

def _fig01b_one_conformer(work, conf):
    g = work[work["conformer"].eq(conf)].copy()
    table = g.pivot_table(index="type_simple", columns="substituent", values="E2_table", aggfunc="sum", fill_value=0)
    table = table.reindex(index=EXPECTED_MAJOR_TYPES, columns=SUB_ORDER).fillna(0)
    return table.reset_index().rename(columns={"type_simple": "Interaction Type"})

def _fig02_long(work):
    idx = pd.MultiIndex.from_product([CONFORMER_ORDER, EXPECTED_MAJOR_TYPES], names=["Conformer", "Interaction Type"])
    out = work.groupby(["conformer", "type_simple"], observed=False).agg(**{"Total E2 (kcal/mol)": ("E2_table", "sum")}).reset_index()
    out = out.rename(columns={"conformer": "Conformer", "type_simple": "Interaction Type"})
    out = out.set_index(["Conformer", "Interaction Type"]).reindex(idx).fillna(0).reset_index()
    return out

def _clean_display_df(df_in):
    out = tidy_round(df_in.copy() if isinstance(df_in, pd.DataFrame) else pd.DataFrame())
    return out if not out.empty else pd.DataFrame({"Note": ["No data available in this run."]})

def _write_title(ws, fmts, row, title, ncols):
    ws.merge_range(row, 0, row, max(ncols - 1, 0), title, fmts["title"])
    ws.set_row(row, 24)

def _write_note(ws, fmts, row, note, ncols):
    if note:
        ws.merge_range(row, 0, row, max(ncols - 1, 0), note, fmts["note"])
        ws.set_row(row, 30)
        return row + 1
    return row

def _write_df(ws, fmts, start_row, start_col, title, df_in):
    df_out = _clean_display_df(df_in)
    ncols = max(len(df_out.columns), 1)
    ws.merge_range(start_row, start_col, start_row, start_col + ncols - 1, title, fmts["section"])
    header_row = start_row + 1
    for j, col in enumerate(df_out.columns):
        ws.write(header_row, start_col + j, str(col), fmts["header"])
        col_vals = df_out[col].astype(str) if col in df_out else pd.Series([], dtype=str)
        width = min(max(11, len(str(col)) + 2, int(col_vals.str.len().quantile(0.90)) + 2 if len(col_vals) else 11), 34)
        ws.set_column(start_col + j, start_col + j, width)
    for i, (_, r) in enumerate(df_out.iterrows(), start=header_row + 1):
        for j, col in enumerate(df_out.columns):
            val = r[col]
            cell_fmt = fmts["num"] if isinstance(val, (int, float, np.integer, np.floating)) and not isinstance(val, bool) and pd.notna(val) else fmts["text"]
            if pd.isna(val):
                ws.write_blank(i, start_col + j, None, fmts["blank"])
            elif cell_fmt is fmts["num"]:
                ws.write_number(i, start_col + j, float(val), cell_fmt)
            else:
                ws.write(i, start_col + j, str(val), cell_fmt)
    ws.autofilter(header_row, start_col, header_row + len(df_out), start_col + ncols - 1)
    return header_row + len(df_out) + 3

def _write_side_by_side(ws, fmts, row, left_title, left_df, right_title, right_df):
    left = _clean_display_df(left_df)
    right = _clean_display_df(right_df)
    left_cols = max(len(left.columns), 1)
    _write_df(ws, fmts, row, 0, left_title, left)
    _write_df(ws, fmts, row, left_cols + 2, right_title, right)
    return row + max(len(left), len(right)) + 4

def _write_fig01_sheet(writer, used, all_df, core_df):
    wb = writer.book
    ws = wb.add_worksheet(_safe_sheet_name("Fig01_E2", used))
    writer.sheets[ws.name] = ws
    fmts = _paper_formats(wb)
    _write_title(ws, fmts, 0, "Figure 1 table: total E2 by interaction type/family/substituent", 21)
    row = _write_note(ws, fmts, 1, "Total E2 only. All interactions use every interaction in mainframe with E2=0 when absent; core interactions restrict atom labels to ≤27.", 21) + 1
    for typ in EXPECTED_MAJOR_TYPES:
        row = _write_side_by_side(
            ws, fmts, row,
            f"All interactions: {typ} total E2", _family_matrix(all_df, typ),
            f"Core interactions: {typ} total E2", _family_matrix(core_df, typ)
        )
    ws.freeze_panes(3, 0)

def _write_fig01b_sheet(writer, used, all_df, core_df):
    wb = writer.book
    ws = wb.add_worksheet(_safe_sheet_name("Fig01_B", used))
    writer.sheets[ws.name] = ws
    fmts = _paper_formats(wb)
    _write_title(ws, fmts, 0, "Figure 1B table: interaction-type total E2 by conformer and substituent", 21)
    row = _write_note(ws, fmts, 1, "Eight compact tables: four conformers for all interactions and the matching four conformers for core interactions only (≤27).", 21) + 1
    for conf in CONFORMER_ORDER:
        row = _write_side_by_side(
            ws, fmts, row,
            f"All interactions: {conf} total E2", _fig01b_one_conformer(all_df, conf),
            f"Core interactions: {conf} total E2", _fig01b_one_conformer(core_df, conf)
        )
    ws.freeze_panes(3, 0)

def _write_fig02_sheet(writer, used, all_df, core_df):
    wb = writer.book
    ws = wb.add_worksheet(_safe_sheet_name("Fig02", used))
    writer.sheets[ws.name] = ws
    fmts = _paper_formats(wb)
    _write_title(ws, fmts, 0, "Figure 2 table: total E2 by interaction type and conformer", 8)
    row = _write_note(ws, fmts, 1, "Total E2 only; absent interactions contribute 0.", 8) + 1
    _write_side_by_side(ws, fmts, row, "All interactions", _fig02_long(all_df), f"Core interactions only (≤{CORE_MAX_ATOM})", _fig02_long(core_df))
    ws.freeze_panes(3, 0)

def _write_partition_sheet(writer, used, sheet_name, partition, description):
    wb = writer.book
    ws = wb.add_worksheet(_safe_sheet_name(sheet_name, used))
    writer.sheets[ws.name] = ws
    fmts = _paper_formats(wb)
    all_df = _scope_df(core_only=False, partition=partition)
    core_df = _scope_df(core_only=True, partition=partition)
    _write_title(ws, fmts, 0, sheet_name, 21)
    row = _write_note(ws, fmts, 1, description, 21) + 1
    for typ in EXPECTED_MAJOR_TYPES:
        row = _write_side_by_side(
            ws, fmts, row,
            f"All interactions: {typ} total E2", _conformer_matrix(all_df, typ),
            f"Core interactions: {typ} total E2", _conformer_matrix(core_df, typ)
        )
    ws.freeze_panes(3, 0)


def _p1p2_metric_table(work, substituent=None, mode="single"):
    g = work[work["conformer"].isin(["P1", "P2"])].copy()
    if substituent is not None:
        g = g[g["substituent"].eq(substituent)].copy()
    if g.empty:
        return pd.DataFrame(columns=["Interaction Type", "P1 total E2", "P2 total E2", "Δ total E2", "P1 count", "P2 count", "Δ count", "P1 mean E2", "P2 mean E2", "Δ mean E2"])

    by_cols = ["substituent", "conformer", "type_simple"] if mode == "mean" else ["conformer", "type_simple"]
    agg = g.groupby(by_cols, observed=False).agg(total_E2=("E2_table", "sum"), count=("count_table", "sum")).reset_index()
    agg["mean_E2"] = agg["total_E2"] / agg["count"].replace(0, np.nan)
    if mode == "mean":
        agg = agg.groupby(["conformer", "type_simple"], observed=False).agg(total_E2=("total_E2", "mean"), count=("count", "mean"), mean_E2=("mean_E2", "mean")).reset_index()

    idx = pd.Index(EXPECTED_MAJOR_TYPES, name="Interaction Type")
    pieces = []
    for metric, label in [("total_E2", "total E2"), ("count", "count"), ("mean_E2", "mean E2")]:
        piv = agg.pivot_table(index="type_simple", columns="conformer", values=metric, aggfunc="sum" if mode != "mean" else "mean", fill_value=0)
        piv = piv.reindex(index=EXPECTED_MAJOR_TYPES, columns=["P1", "P2"]).fillna(0)
        out = pd.DataFrame(index=idx)
        out[f"P1 {label}"] = piv["P1"].values
        out[f"P2 {label}"] = piv["P2"].values
        out[f"Δ {label}"] = piv["P1"].values - piv["P2"].values
        pieces.append(out)
    return pd.concat(pieces, axis=1).reset_index()


def _write_p1p2_sheet(writer, used):
    wb = writer.book
    ws = wb.add_worksheet(_safe_sheet_name("P1_vs_P2", used))
    writer.sheets[ws.name] = ws
    fmts = _paper_formats(wb)
    all_df = _scope_df(core_only=False)
    _write_title(ws, fmts, 0, "P1 vs P2 by interaction type", 10)
    note = (
        "Eleven tables on one sheet. Nine substituent-specific tables plus one total table and one mean-across-substituents table. "
        "Metrics are total E2, count of existing interactions, and mean E2 of existing interactions. "
        "All interactions are used here (full 73-interaction universe with absent interactions contributing zero to totals and counts)."
    )
    row = _write_note(ws, fmts, 1, note, 10) + 1
    for sub in SUB_ORDER:
        row = _write_df(ws, fmts, row, 0, f"{sub}: P1 vs P2", _p1p2_metric_table(all_df, substituent=sub, mode="single"))
    row = _write_df(ws, fmts, row, 0, "Total across all substituents: P1 vs P2", _p1p2_metric_table(all_df, substituent=None, mode="single"))
    _write_df(ws, fmts, row, 0, "Mean across substituents: P1 vs P2", _p1p2_metric_table(all_df, substituent=None, mode="mean"))
    ws.freeze_panes(3, 0)

def _safe_write_paper_workbook(output_path):
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    try:
        export_paper_tables_workbook(output_path)
        return output_path
    except Exception as exc:
        print(f"Paper tables workbook skipped due to error: {exc}")
        return None

def export_paper_tables_workbook(output_path):
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    used = set()
    all_df = _scope_df(core_only=False)
    core_df = _scope_df(core_only=True)
    with pd.ExcelWriter(output_path, engine="xlsxwriter") as writer:
        _write_fig01_sheet(writer, used, all_df, core_df)
        _write_fig01b_sheet(writer, used, all_df, core_df)
        _write_fig02_sheet(writer, used, all_df, core_df)
        _write_partition_sheet(
            writer, used, "NMe2 - Propanol Specific", "nonring",
            "Excluded pure 1–9 internal labels such as C1 or C1/C2. Outside connections such as C1/C19 are retained."
        )
        _write_partition_sheet(
            writer, used, "Ring Specific", "ring",
            "Kept only interactions where donor or acceptor labels are pure 1–9 internal labels such as C1 or C1/C2."
        )
        _write_p1p2_sheet(writer, used)
    print(f"Paper tables workbook written: {output_path}")
    return output_path

_safe_write_paper_workbook(os.path.join(PAPERTABLEDIR, PAPER_TABLE_WORKBOOK))

type_conf = type_sample.groupby(["conformer", "conformer_family", "type_simple"], observed=False).agg(total_E2=("total_E2", "sum"), count_existing=("count_existing", "sum"), mean_E2=("mean_E2", "mean"), max_E2=("max_E2", "max")).reset_index()
type_conf["fraction_of_conformer_E2"] = type_conf["total_E2"] / type_conf.groupby("conformer")["total_E2"].transform("sum").replace(0, np.nan)


pn_total = summary.groupby(["substituent", "sigma_p"], observed=False).apply(lambda g: pd.Series({
    "P_mean_total_E2": g.loc[g["conformer_family"].eq("P"), "total_E2"].mean(),
    "N_mean_total_E2": g.loc[g["conformer_family"].eq("N"), "total_E2"].mean(),
    "P_mean_npi_E2": g.loc[g["conformer_family"].eq("P"), "npi_total_E2"].mean(),
    "N_mean_npi_E2": g.loc[g["conformer_family"].eq("N"), "npi_total_E2"].mean(),
    "P_mean_non_npi_E2": g.loc[g["conformer_family"].eq("P"), "non_npi_total_E2"].mean(),
    "N_mean_non_npi_E2": g.loc[g["conformer_family"].eq("N"), "non_npi_total_E2"].mean(),
})).reset_index()
pn_total["P_minus_N_total_E2"] = pn_total["P_mean_total_E2"] - pn_total["N_mean_total_E2"]
pn_total["P_minus_N_npi_E2"] = pn_total["P_mean_npi_E2"] - pn_total["N_mean_npi_E2"]
pn_total["P_minus_N_non_npi_E2"] = pn_total["P_mean_non_npi_E2"] - pn_total["N_mean_non_npi_E2"]

idx = pd.MultiIndex.from_product([SUB_ORDER, CONFORMER_ORDER], names=["substituent", "conformer"])
E = df.pivot_table(index=["substituent", "conformer"], columns="interaction_id", values="E2_eff", aggfunc="mean").reindex(idx)

core_interaction_ids = dictionary.loc[dictionary["core_interaction"].eq(True), "interaction_id"].tolist()
E_core = E[[c for c in E.columns if c in core_interaction_ids]].copy()

geom_vars = [c for c in ["distance", "angle", "dihedral"] if c in df.columns]
existing = df[df["possibility_flag"].eq(1) & df["exists_flag"].eq(1) & df["E2_present"].notna()].copy()
geom_rows = []
for typ in type_order(existing["type_simple"]):
    for fam in ["all", "N", "P"]:
        g = existing[existing["type_simple"].eq(typ)] if fam == "all" else existing[existing["type_simple"].eq(typ) & existing["conformer_family"].eq(fam)]
        for var in geom_vars:
            r, n = corr_pair(g[var], g["E2_present"])
            geom_rows.append({"type_simple": typ, "conformer_family": fam, "geometry_variable": var, "n_points": n, "corr_E2_vs_geometry": r, "interpretation_hint": "E2 rises with coordinate" if pd.notna(r) and r > 0.35 else "E2 rises as coordinate decreases" if pd.notna(r) and r < -0.35 else "weak/mixed"})
geom_corr = pd.DataFrame(geom_rows).sort_values("corr_E2_vs_geometry", key=lambda s: s.abs(), ascending=False)


pareto = df[df["possibility_flag"].eq(1)].groupby("interaction_id").agg(interaction_code=("interaction_code", "first"), type_simple=("type_simple", "first"), total_E2=("E2_eff", "sum"), max_E2=("E2_present", "max"), n_exists=("exists_flag", "sum")).reset_index().sort_values("total_E2", ascending=False)
pareto["long_label"] = pareto["interaction_id"].map(label_map)
pareto["rank"] = np.arange(1, len(pareto) + 1)
pareto["cumulative_E2"] = pareto["total_E2"].cumsum()
pareto["cumulative_fraction"] = pareto["cumulative_E2"] / pareto["total_E2"].sum()

ordered_sample_index = pd.MultiIndex.from_product([SUB_ORDER, CONFORMER_ORDER], names=["substituent", "conformer"])

def make_type_linegrid(source_df, value_col, outname, main_title, family_col="conformer_family", families=("N", "P"), type_filter=None, y_label="ΣE2"):
    work = source_df.copy()
    if type_filter is not None:
        work = work[work["type_simple"].isin(type_filter)].copy()
    tlist = type_order(work["type_simple"])
    if not tlist:
        return
    ncols = 2
    nrows = math.ceil(len(tlist) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(14, max(4 * nrows, 6)), sharex=True)
    axes = np.array(axes).reshape(-1)
    for ax, typ in zip(axes, tlist):
        g = work[work["type_simple"].eq(typ)].copy()
        for fam in families:
            gg = g[g[family_col].eq(fam)].sort_values("sigma_p")
            if not gg.empty:
                ax.plot(gg["substituent"], gg[value_col], marker="o", label=fam)
        ax.set_title(typ)
        ax.set_ylabel(y_label)
        ax.grid(alpha=0.25)
    for ax in axes[len(tlist):]:
        ax.axis("off")
    for ax in axes[:len(tlist)]:
        ax.set_xlabel("Substituent, EDG → EWG")
    handles, labels = axes[0].get_legend_handles_labels()
    if handles:
        axes[0].legend(handles, labels, title="Family", fontsize=8)
    fig.suptitle(main_title, y=1.02, fontsize=14)
    savefig(outname)

fig01_all = type_sample.groupby(["substituent", "sigma_p", "conformer_family", "type_simple"], observed=False).agg(
    mean_total_E2=("total_E2", "mean"),
    mean_count_existing=("count_existing", "mean")
).reset_index()
fig01_core = type_sample_core.groupby(["substituent", "sigma_p", "conformer_family", "type_simple"], observed=False).agg(
    mean_total_E2=("total_E2", "mean"),
    mean_count_existing=("count_existing", "mean")
).reset_index()

make_type_linegrid(fig01_all.rename(columns={"mean_total_E2": "value"}), "value", "FIG_01_reference_style_comparable_energies_N_vs_P.png", "Comparable interaction energies by interaction type (all atoms)", y_label="Mean ΣE2")
make_type_linegrid(fig01_all.rename(columns={"mean_count_existing": "value"}), "value", "FIG_01A_reference_style_interaction_counts_N_vs_P.png", "Mean number of interactions by interaction type (all atoms)", y_label="Mean count")
make_type_linegrid(fig01_core.rename(columns={"mean_total_E2": "value"}), "value", "FIG_01_core_atoms_comparable_energies_N_vs_P.png", f"Comparable interaction energies by interaction type (core atoms only, atom index ≤ {CORE_MAX_ATOM})", y_label="Mean ΣE2")
make_type_linegrid(fig01_core.rename(columns={"mean_count_existing": "value"}), "value", "FIG_01A_core_atoms_interaction_counts_N_vs_P.png", f"Mean number of interactions by interaction type (core atoms only, atom index ≤ {CORE_MAX_ATOM})", y_label="Mean count")

def plot_fig01b_variant(source_df, filename, title_suffix="", exclude_types=None, zoom=False):
    exclude_types = set(exclude_types or [])
    fig, axes = plt.subplots(2, 2, figsize=(15, 9), sharex=True, sharey=False)
    axes = np.array(axes).reshape(-1)
    all_plotted_types = []
    for ax, conf in zip(axes, CONFORMER_ORDER):
        g_conf = source_df[source_df["conformer"].eq(conf)].copy()
        g_conf = g_conf[~g_conf["type_simple"].isin(exclude_types)].copy()
        plotted_vals = []
        for typ in type_order(g_conf["type_simple"]):
            g = g_conf[g_conf["type_simple"].eq(typ)].sort_values("sigma_p")
            if not g.empty:
                ax.plot(g["substituent"], g["total_E2"], marker="o", label=typ)
                plotted_vals.extend(g["total_E2"].dropna().tolist())
                all_plotted_types.append(typ)
        if zoom and plotted_vals:
            lo, hi = min(plotted_vals), max(plotted_vals)
            pad = max(0.06 * (hi - lo), 0.05) if hi > lo else max(0.08 * max(abs(hi), 1.0), 0.05)
            ax.set_ylim(lo - pad, hi + pad)
        ax.set_title(conf)
        ax.set_xlabel("Substituent, EDG → EWG")
        ax.set_ylabel("ΣE2")
        ax.grid(alpha=0.25)
    handles, labels = axes[0].get_legend_handles_labels()
    if handles:
        seen = set()
        h2, l2 = [], []
        for h, l in zip(handles, labels):
            if l not in seen:
                seen.add(l)
                h2.append(h)
                l2.append(l)
        axes[0].legend(h2, l2, title="Interaction type", fontsize=8)
    suffix = f" ({title_suffix})" if title_suffix else ""
    fig.suptitle(f"E2 by substituent with interaction-type lines{suffix}", y=1.02, fontsize=14)
    savefig(filename)

plot_fig01b_variant(
    type_sample,
    "FIG_01B_E2_by_substituent_type_lines_each_conformer.png"
)
plot_fig01b_variant(
    type_sample,
    "FIG_01B2_E2_by_substituent_type_lines_each_conformer_no_nsig_ssig.png",
    title_suffix="without n→σ* and σ→σ*",
    exclude_types={"n→σ*", "σ→σ*"},
    zoom=True
)
plot_fig01b_variant(
    type_sample,
    "FIG_01B3_E2_by_substituent_type_lines_each_conformer_no_npi.png",
    title_suffix="without n→π*",
    exclude_types={"n→π*"},
    zoom=True
)
plot_fig01b_variant(
    type_sample_core,
    "FIG_01B4_E2_by_substituent_type_lines_each_conformer_core_atoms.png",
    title_suffix="core atoms only (≤27)"
)
plot_fig01b_variant(
    type_sample_core,
    "FIG_01B5_E2_by_substituent_type_lines_each_conformer_core_atoms_no_nsig_ssig.png",
    title_suffix="core atoms only (≤27), without n→σ* and σ→σ*",
    exclude_types={"n→σ*", "σ→σ*"},
    zoom=True
)
plot_fig01b_variant(
    type_sample_core,
    "FIG_01B6_E2_by_substituent_type_lines_each_conformer_core_atoms_no_npi.png",
    title_suffix="core atoms only (≤27), without n→π*",
    exclude_types={"n→π*"},
    zoom=True
)

def plot_fig02_variant(source_df, filename, title_suffix=""):
    type_conf_local = source_df.groupby(
        ["conformer", "conformer_family", "type_simple"], observed=False
    ).agg(
        total_E2=("total_E2", "sum"),
        mean_E2=("mean_E2", "mean"),
        max_E2=("max_E2", "max")
    ).reset_index()
    if type_conf_local.empty:
        return
    fig, ax = plt.subplots(figsize=(10, 5))
    tlist = type_order(type_conf_local["type_simple"])
    for conf in CONFORMER_ORDER:
        g = type_conf_local[type_conf_local["conformer"].eq(conf)].set_index("type_simple").reindex(tlist)
        ax.plot(tlist, g["total_E2"], marker="o", label=conf)
    ttl = "Interaction-type stabilization by conformer"
    if title_suffix:
        ttl += f" ({title_suffix})"
    ax.set_title(ttl)
    ax.set_xlabel("Interaction type")
    ax.set_ylabel("ΣE2")
    ax.grid(alpha=0.25)
    ax.legend(title="Conformer")
    plt.xticks(rotation=25, ha="right")
    savefig(filename)

plot_fig02_variant(type_sample, "FIG_02_interaction_type_by_conformer.png")
plot_fig02_variant(type_sample_core, "FIG_02A_interaction_type_by_conformer_core_atoms.png", "core atoms only, ≤27")


non_types = [t for t in type_order(type_sample["type_simple"]) if t != "n→π*"]
if non_types:

    def _type_slug(typ):
        return (
            sstr(typ).replace("→", "_to_").replace("*", "star")
            .replace("π", "pi").replace("σ", "sigma")
            .replace("/", "_").replace(" ", "_").replace("-", "_")
        )

    def _fig03a_presentation(src, scope_label, typ, filename):
        g = src[src["type_simple"].eq(typ)].copy()
        if g.empty:
            return None
        x = np.arange(len(SUB_ORDER))
        fig, ax = plt.subplots(figsize=(12, 7))
        plotted = []
        for conf in CONFORMER_ORDER:
            gg = g[g["conformer"].eq(conf)].set_index("substituent").reindex(SUB_ORDER)
            y = gg["total_E2"].astype(float)
            if y.notna().any():
                ax.plot(x, y, marker="o", linewidth=2.2, markersize=6, label=conf)
                plotted.extend(y.dropna().tolist())
                last_valid = y.dropna()
                if not last_valid.empty:
                    last_sub = last_valid.index[-1]
                    last_x = SUB_ORDER.index(last_sub)
                    last_y = float(last_valid.iloc[-1])
                    ax.annotate(f"{conf}", (last_x, last_y), xytext=(6, 0), textcoords="offset points", fontsize=8, va="center")
        if plotted:
            lo, hi = min(plotted), max(plotted)
            pad = max((hi - lo) * 0.10, 0.05) if hi > lo else max(abs(hi) * 0.10, 0.05)
            ax.set_ylim(lo - pad, hi + pad)
        ax.set_xticks(x)
        ax.set_xticklabels(SUB_ORDER, rotation=35, ha="right")
        ax.set_title(f"FIG 03A presentation — {scope_label}: {typ}", fontsize=14, pad=14)
        ax.set_xlabel("Substituent, EDG → EWG")
        ax.set_ylabel(f"ΣE2 for {typ}")
        ax.grid(axis="y", alpha=0.28)
        ax.legend(title="Conformer", ncol=4, loc="upper center", bbox_to_anchor=(0.5, -0.18), frameon=False)
        add_hammett_top_axis(ax, SUB_ORDER)
        return savefig_presentation(filename)

    for _scope_i, (_src, _scope_name, _scope_slug) in enumerate([
        (type_sample, "All interactions", "all_interactions"),
        (type_sample_core, f"Core interactions only (≤{CORE_MAX_ATOM})", "core_interactions"),
    ], start=1):
        for _type_i, _typ in enumerate(non_types, start=1):
            _fig03a_presentation(
                _src, _scope_name, _typ,
                f"FIG_03A_P{_scope_i}_{_type_i}_{_scope_slug}_{_type_slug(_typ)}.png"
            )

    def _line_detail_graph(src_df, conf, typ, filename, scope_label):
        g = src_df[(src_df["conformer"].eq(conf)) & (src_df["type_simple"].eq(typ))].copy()
        if g.empty:
            return None
        g = g.set_index("substituent").reindex(SUB_ORDER).reset_index()
        x = np.arange(len(SUB_ORDER), dtype=float)
        y = pd.to_numeric(g["total_E2"], errors="coerce").astype(float)
        fig, ax = plt.subplots(figsize=(10.5, 6.2))
        ax.plot(x, y, marker="o", linewidth=2.5, markersize=7)
        for xi, yi in zip(x, y):
            if pd.notna(yi):
                ax.annotate(f"{yi:.2f}", (xi, yi), xytext=(0, 7), textcoords="offset points", ha="center", fontsize=8)

        valid = pd.DataFrame({"x_plot": x, "y": y}).dropna()
        if len(valid) >= 3:
            r = float(np.corrcoef(valid["x_plot"], valid["y"])[0, 1]) if valid["x_plot"].nunique() > 1 and valid["y"].nunique() > 1 else np.nan
            if pd.notna(r) and (r ** 2) > LINE_DETAIL_R_THRESHOLD:
                slope, intercept = np.polyfit(valid["x_plot"], valid["y"], 1)
                x_fit = np.array([valid["x_plot"].min(), valid["x_plot"].max()])
                y_fit = slope * x_fit + intercept
                ax.plot(x_fit, y_fit, linestyle="--", linewidth=1.7, alpha=0.90)
                sign = "+" if slope >= 0 else "−"
                ax.text(
                    0.02, 0.98,
                    f"Strong EDG→EWG trend: R² = {r**2:.2f} ({sign})",
                    transform=ax.transAxes, va="top", ha="left", fontsize=9,
                    bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "alpha": 0.86, "edgecolor": "#BBBBBB"}
                )

        ax.set_xticks(x)
        ax.set_xticklabels(SUB_ORDER, rotation=35, ha="right")
        ax.set_xlabel("Substituent, EDG → EWG")
        ax.set_ylabel("Total E2 (kcal/mol)")
        ax.set_title(f"{conf} — {typ} across the substituent series ({scope_label})")
        ax.grid(axis="y", alpha=0.30)
        add_hammett_top_axis(ax, SUB_ORDER)
        return savefig_presentation(filename)

    for _conf in CONFORMER_ORDER:
        for _typ in EXPECTED_MAJOR_TYPES:
            _line_detail_graph(type_sample, _conf, _typ, f"FIG_01B20_{_conf}_{_type_slug(_typ)}.png", "all interactions")
            _line_detail_graph(type_sample_core, _conf, _typ, f"FIG_01B20_core_{_conf}_{_type_slug(_typ)}.png", f"core interactions only (≤{CORE_MAX_ATOM})")

    def _slide_type_label(typ):
        return sstr(typ).replace("→", "-")

    # kept seperate so line styles still work in b/w
    PAPER_PANEL_COLORS = {"N1": "#1f77b4", "N2": "#ff7f0e", "P1": "#2ca02c", "P2": "#d62728"}
    PAPER_LINESTYLES = {"N1": "-", "N2": "--", "P1": "-.", "P2": ":"}
    PAPER_MARKERS = {"N1": "o", "N2": "s", "P1": "^", "P2": "D"}

    def _style_paper_axes(ax):
        ax.tick_params(axis="both", which="both", direction="in", width=0.75, labelsize=8, top=True, right=True)
        for spine in ax.spines.values():
            spine.set_linewidth(1.5)

    def _trend_panel(ax, src, typ, lgnd):
        g = src[src["type_simple"].eq(typ)].copy()
        x = np.arange(len(SUB_ORDER), dtype=float)
        yy = []
        for conf in CONFORMER_ORDER:
            gg = g[g["conformer"].eq(conf)].set_index("substituent").reindex(SUB_ORDER)
            y = pd.to_numeric(gg["total_E2"], errors="coerce").astype(float)
            ok = y.dropna()
            if ok.empty or np.allclose(ok.to_numpy(dtype=float), 0.0):
                continue
            ax.plot(
                x, y,
                color=PAPER_PANEL_COLORS.get(conf),
                linestyle=PAPER_LINESTYLES.get(conf, "-"),
                marker=PAPER_MARKERS.get(conf, "o"),
                markersize=6.6,
                linewidth=2.1,
                label=conf,
                markeredgecolor="white",
                markeredgewidth=0.75,
            )
            yy.extend(ok.tolist())

        if yy:
            lo, hi = min(yy), max(yy)
            pad = max(0.12 * (hi - lo), 0.04) if hi > lo else max(abs(hi) * 0.14, 0.05)
            ax.set_ylim(lo - pad, hi + pad)
        ax.set_xticks(x)
        ax.set_xticklabels(SUB_ORDER, rotation=32, ha="right", fontsize=8)
        ax.grid(axis="y", alpha=0.24, linewidth=0.8)
        ax.set_xlabel("Substituent", fontsize=9)
        ax.set_ylabel("Total E2 (kcal/mol)", fontsize=9)
        _style_paper_axes(ax)
        handles, labels = ax.get_legend_handles_labels()
        if handles:
            ax.legend(handles, labels, loc="best", frameon=True, framealpha=0.95, fontsize=8, title=lgnd, title_fontsize=8, borderpad=0.35, handlelength=2.2)

    def _trend_figure():
        panels = [
            ("n→σ*", "n-σ*"),
            ("π→σ*", "π-σ*"),
            ("σ→σ*", "σ-σ*"),
            ("σ→π*", "σ-π*"),
        ]
        fig, axes = plt.subplots(2, 2, figsize=(16.4, 11.4), facecolor="white")
        axes = np.array(axes).reshape(-1)
        for ax, (typ, label) in zip(axes, panels):
            _trend_panel(ax, type_sample_core, typ, label)
        plt.subplots_adjust(left=0.07, right=0.985, top=0.985, bottom=0.10, wspace=0.18, hspace=0.28)
        savefig_paper("FIG_01B_core_paper_4panel.png")

    _trend_figure()

    GEOM_LABELS = {"distance": "Distance (Å)", "angle": "Angle (°)", "dihedral": "Dihedral (°)"}
    CONFORMER_COLORS = {"N1": "#1f77b4", "N2": "#ff7f0e", "P1": "#2ca02c", "P2": "#d62728"}
    GEOM_BUBBLE_CASES = ["n→π*", "n→σ* CH", "n→σ* CC", "π→σ*", "σ→π*", "σ→σ*"]

    def _paper_geometry_case(row):
        typ = row.get("type_simple", "")
        if typ == "n→π*":
            return "n→π*"
        if typ == "n→σ*":
            acc_atoms = atom_tokens(row.get("acceptor_labels", ""))
            return "n→σ* CH" if any(atom_element(a) == "H" for a in acc_atoms) else "n→σ* CC"
        if typ in {"π→σ*", "σ→π*", "σ→σ*"}:
            return typ
        return None

    def _geom_case_slug(case):
        return {
            "n→π*": "n_to_pistar",
            "n→σ* CH": "n_to_sigmastar_CH",
            "n→σ* CC": "n_to_sigmastar_CC",
            "π→σ*": "pi_to_sigmastar",
            "σ→π*": "sigma_to_pistar",
            "σ→σ*": "sigma_to_sigmastar",
        }.get(case, _type_slug(case))

    def _legend_case_label(case):
        return case.replace("→", "-")

    def _sigma_pi_ch_to_co(row):
        donor_elems = [atom_element(a) for a in atom_tokens(row.get("donor_labels", ""))]
        acceptor_elems = [atom_element(a) for a in atom_tokens(row.get("acceptor_labels", ""))]
        return ("C" in donor_elems and "H" in donor_elems and "C" in acceptor_elems and "O" in acceptor_elems)

    def _filter_case_rows(source_df, case):
        work = source_df.copy()
        work["_paper_geometry_case"] = work.apply(_paper_geometry_case, axis=1)
        work = work[work["_paper_geometry_case"].eq(case)].copy()
        if case == "σ→π*":
            work = work[work.apply(_sigma_pi_ch_to_co, axis=1)].copy()
        return work

    def _geometry_pair_for_case(g, case):
        preferred = {
            "n→π*": [("distance", "angle"), ("distance", "E2_present"), ("angle", "E2_present")],
            "n→σ* CH": [("distance", "dihedral"), ("distance", "angle"), ("distance", "E2_present")],
            "n→σ* CC": [("distance", "angle"), ("angle", "distance"), ("angle", "E2_present")],
            "π→σ*": [("distance", "angle"), ("distance", "dihedral"), ("distance", "E2_present"), ("angle", "E2_present")],
            "σ→π*": [("distance", "dihedral")],
            "σ→σ*": [("distance", "angle"), ("distance", "dihedral"), ("angle", "E2_present"), ("distance", "E2_present")],
        }.get(case, [("distance", "angle"), ("distance", "dihedral"), ("angle", "dihedral"), ("angle", "E2_present")])
        for xvar, yvar in preferred:
            if xvar in g.columns and yvar in g.columns:
                need_cols = [xvar, yvar] + (["E2_present"] if "E2_present" not in {xvar, yvar} else [])
                n = g[need_cols].replace([np.inf, -np.inf], np.nan).dropna().shape[0]
                if n > 0:
                    return xvar, yvar
        return None

    NPI_WINDOW = {"distance_max": 3.2, "angle_min": 95.0, "angle_max": 125.0}
    NSIG_CH_WINDOW = {"distance_min": 2.2, "distance_max": 2.6, "abs_dihedral_max": 53.0}
    NSIG_CC_WINDOW = {"angle_min": 130.0, "angle_max": 180.0}

    def _abs_wrap_dihedral(val):
        if pd.isna(val):
            return np.nan
        x = abs(float(val))
        return min(x, abs(360.0 - x)) if x > 180 else x

    def _within_predefined_window(row, case):
        dist = pd.to_numeric(pd.Series([row.get("distance")]), errors="coerce").iloc[0]
        ang = pd.to_numeric(pd.Series([row.get("angle")]), errors="coerce").iloc[0]
        dih = pd.to_numeric(pd.Series([row.get("dihedral")]), errors="coerce").iloc[0]
        if case == "n→π*":
            return pd.notna(dist) and pd.notna(ang) and dist <= NPI_WINDOW["distance_max"] and NPI_WINDOW["angle_min"] <= ang <= NPI_WINDOW["angle_max"]
        if case == "n→σ* CH":
            ad = _abs_wrap_dihedral(dih)
            return pd.notna(dist) and pd.notna(ad) and NSIG_CH_WINDOW["distance_min"] <= dist <= NSIG_CH_WINDOW["distance_max"] and ad <= NSIG_CH_WINDOW["abs_dihedral_max"]
        if case == "n→σ* CC":
            return pd.notna(ang) and NSIG_CC_WINDOW["angle_min"] <= ang <= NSIG_CC_WINDOW["angle_max"]
        return False

    def _within_dynamic_window(row, active_df, pair):
        geom_vars = [v for v in pair if v in {"distance", "angle", "dihedral"}]
        if not geom_vars:
            return True
        ref = active_df.replace([np.inf, -np.inf], np.nan).dropna(subset=geom_vars)
        if ref.empty:
            return False
        for var in geom_vars:
            val = pd.to_numeric(pd.Series([row.get(var)]), errors="coerce").iloc[0]
            if pd.isna(val):
                return False
            lo = float(ref[var].min())
            hi = float(ref[var].max())
            span = hi - lo
            pad = max(0.05 * span, 1e-6)
            if val < lo - pad or val > hi + pad:
                return False
        return True

    def _filter_window_fit_absent(g_zero, g_active, case, pair):
        if g_zero.empty:
            return g_zero.copy()
        if case in {"n→π*", "n→σ* CH", "n→σ* CC"}:
            return g_zero[g_zero.apply(lambda r: _within_predefined_window(r, case), axis=1)].copy()
        return g_zero[g_zero.apply(lambda r: _within_dynamic_window(r, g_active, pair), axis=1)].copy()

    def _bubble_size(vals, vmin=None, vmax=None):
        v = pd.to_numeric(vals, errors="coerce").fillna(0).clip(lower=0)
        if len(v) == 0:
            return np.array([], dtype=float)
        lo = float(v.min()) if vmin is None or pd.isna(vmin) else float(vmin)
        hi = float(v.max()) if vmax is None or pd.isna(vmax) else float(vmax)
        if hi <= lo:
            return np.repeat(180.0, len(v))
        return 70 + 560 * (v - lo) / (hi - lo)

    def _limit_from_series(values, include_zero=False):
        s = pd.to_numeric(pd.Series(values), errors="coerce").replace([np.inf, -np.inf], np.nan).dropna()
        if s.empty:
            return None
        lo = float(s.min())
        hi = float(s.max())
        if include_zero:
            lo = min(lo, 0.0)
            hi = max(hi, 0.0)
        if math.isclose(lo, hi):
            pad = max(abs(lo) * 0.10, 1.0)
        else:
            pad = max((hi - lo) * 0.08, 1e-6)
        return (lo - pad, hi + pad)

    def _mk_bins(values, geom_var):
        vals = pd.to_numeric(pd.Series(values), errors="coerce").replace([np.inf, -np.inf], np.nan).dropna()
        if vals.empty:
            return None
        lo = float(vals.min())
        hi = float(vals.max())
        if math.isclose(lo, hi):
            pad = 0.10 if geom_var == "distance" else 5.0
            lo -= pad
            hi += pad
        n_bins = min(12, max(6, int(math.ceil(math.sqrt(len(vals))))))
        return np.linspace(lo, hi, n_bins + 1)

    def _mk_plot_spec(case, source_frames):
        work_all = pd.concat([_filter_case_rows(src, case) for src in source_frames], ignore_index=True, sort=False)
        g_active = work_all[
            work_all["possibility_flag"].eq(1)
            & work_all["exists_flag"].eq(1)
            & work_all["E2_present"].notna()
        ].copy()
        g_zero = work_all[
            work_all["possibility_flag"].eq(1)
            & work_all["exists_flag"].eq(0)
        ].copy()
        pair = _geometry_pair_for_case(g_active, case)
        if pair is None:
            pair = _geometry_pair_for_case(pd.concat([g_active.copy(), g_zero.assign(E2_present=0.0)], ignore_index=True), case)
        if pair is None:
            return {
                "pair": None, "xlim": None, "ylim": None,
                "bubble_min": np.nan, "bubble_max": np.nan,
                "histogram_mode": False, "hist_geom_var": None,
                "hist_bins": None, "hist_xlim": None, "hist_ylim": None,
            }
        xvar, yvar = pair
        histogram_mode = xvar == "E2_present" or yvar == "E2_present"
        geom_var = yvar if xvar == "E2_present" else xvar
        g_zero = _filter_window_fit_absent(g_zero, g_active, case, pair)
        active_clean = g_active.replace([np.inf, -np.inf], np.nan).dropna(subset=[xvar, yvar, "E2_present"])
        if histogram_mode:
            g_zero = g_zero.assign(E2_present=0.0)
        zero_clean = g_zero.replace([np.inf, -np.inf], np.nan).dropna(subset=[xvar, yvar])
        xvals = pd.concat([active_clean[xvar], zero_clean[xvar]], ignore_index=True) if not zero_clean.empty else active_clean[xvar]
        yvals = pd.concat([active_clean[yvar], zero_clean[yvar]], ignore_index=True) if not zero_clean.empty else active_clean[yvar]

        hist_bins = None
        hist_xlim = None
        hist_ylim = None
        if histogram_mode:
            geom_values = pd.concat([
                pd.to_numeric(active_clean[geom_var], errors="coerce"),
                pd.to_numeric(zero_clean[geom_var], errors="coerce"),
            ], ignore_index=True).dropna()
            hist_bins = _mk_bins(geom_values, geom_var)
            if hist_bins is not None:
                hist_xlim = (float(hist_bins[0]), float(hist_bins[-1]))
                max_count = 0
                for src in source_frames:
                    src_case = _filter_case_rows(src, case)
                    src_active = src_case[
                        src_case["possibility_flag"].eq(1)
                        & src_case["exists_flag"].eq(1)
                        & src_case["E2_present"].notna()
                    ].copy()
                    src_zero = src_case[
                        src_case["possibility_flag"].eq(1)
                        & src_case["exists_flag"].eq(0)
                    ].copy()
                    src_zero = _filter_window_fit_absent(src_zero, src_active, case, pair)
                    active_vals = pd.to_numeric(src_active.get(geom_var, pd.Series(dtype=float)), errors="coerce").dropna()
                    zero_vals = pd.to_numeric(src_zero.get(geom_var, pd.Series(dtype=float)), errors="coerce").dropna()
                    if not active_vals.empty:
                        max_count = max(max_count, int(np.histogram(active_vals, bins=hist_bins)[0].max()))
                    if not zero_vals.empty:
                        max_count = max(max_count, int(np.histogram(zero_vals, bins=hist_bins)[0].max()))
                hist_ylim = (0.0, max(1.0, math.ceil(max_count * 1.18)))

        return {
            "pair": pair,
            "xlim": _limit_from_series(xvals, include_zero=(xvar == "E2_present")),
            "ylim": _limit_from_series(yvals, include_zero=(yvar == "E2_present")),
            "bubble_min": float(active_clean["E2_present"].min()) if not active_clean.empty else np.nan,
            "bubble_max": float(active_clean["E2_present"].max()) if not active_clean.empty else np.nan,
            "histogram_mode": histogram_mode,
            "hist_geom_var": geom_var if histogram_mode else None,
            "hist_bins": hist_bins,
            "hist_xlim": hist_xlim,
            "hist_ylim": hist_ylim,
        }

    def _set_lims(ax, spec):
        if spec.get("histogram_mode"):
            if spec.get("hist_xlim") is not None:
                ax.set_xlim(spec["hist_xlim"])
            if spec.get("hist_ylim") is not None:
                ax.set_ylim(spec["hist_ylim"])
            return
        if spec.get("xlim") is not None:
            ax.set_xlim(spec["xlim"])
        if spec.get("ylim") is not None:
            ax.set_ylim(spec["ylim"])

    def _freq_hist(ax, g_active, g_zero, spec, for_panel=False):
        geom_var = spec["hist_geom_var"]
        bins = spec.get("hist_bins")
        if geom_var is None or bins is None:
            return
        centers = (bins[:-1] + bins[1:]) / 2.0
        widths = np.diff(bins) * 0.92
        order = list(CONFORMER_ORDER) + [c for c in pd.unique(g_active["conformer"]) if c not in CONFORMER_ORDER]
        bottom = np.zeros(len(centers), dtype=float)
        plotted_any = False
        for idx, conf in enumerate(order):
            gg = g_active[g_active["conformer"].eq(conf)].copy()
            vals = pd.to_numeric(gg.get(geom_var, pd.Series(dtype=float)), errors="coerce").dropna()
            if vals.empty:
                continue
            counts, _ = np.histogram(vals, bins=bins)
            color = CONFORMER_COLORS.get(conf, f"C{idx % 10}")
            ax.bar(
                centers, counts, width=widths, bottom=bottom,
                color=color, alpha=0.76, edgecolor="white", linewidth=0.45,
                label=conf, align="center",
            )
            bottom += counts
            plotted_any = True

        zero_vals = pd.to_numeric(g_zero.get(geom_var, pd.Series(dtype=float)), errors="coerce").dropna()
        if not zero_vals.empty:
            zero_counts, _ = np.histogram(zero_vals, bins=bins)
            ax.bar(
                centers, zero_counts, width=widths,
                facecolor="none", edgecolor="black", linewidth=0.95,
                hatch="///", label="window-fit, not observed", align="center", zorder=5,
            )
        if not plotted_any and zero_vals.empty:
            ax.text(0.5, 0.5, "No interactions in the selected geometry bins", transform=ax.transAxes, ha="center", va="center", fontsize=9)

    def _geom_or_hist(ax, g_active, g_zero, spec, for_panel=False):
        xvar, yvar = spec["pair"]
        if spec.get("histogram_mode"):
            _freq_hist(ax, g_active, g_zero, spec, for_panel=for_panel)
        else:
            order = list(CONFORMER_ORDER) + [c for c in pd.unique(g_active["conformer"]) if c not in CONFORMER_ORDER]
            for conf in order:
                gg = g_active[g_active["conformer"].eq(conf)].copy()
                if gg.empty:
                    continue
                ax.scatter(
                    gg[xvar], gg[yvar],
                    s=_bubble_size(gg["E2_present"], spec.get("bubble_min"), spec.get("bubble_max")),
                    color=CONFORMER_COLORS.get(conf, None),
                    alpha=0.68,
                    edgecolors="#333333",
                    linewidths=0.35,
                    label=conf,
                )
            if not g_zero.empty:
                ax.scatter(g_zero[xvar], g_zero[yvar], s=18, color="black", alpha=0.95, linewidths=0, label="window-fit, not observed", zorder=5)
        _set_lims(ax, spec)

    def _plot_bits(source_df, case, spec):
        work = _filter_case_rows(source_df, case)
        g_active = work[
            work["possibility_flag"].eq(1)
            & work["exists_flag"].eq(1)
            & work["E2_present"].notna()
        ].copy()
        if spec["pair"] is None:
            g_zero = work.iloc[0:0].copy()
            return g_active, g_zero
        xvar, yvar = spec["pair"]
        g_zero = work[
            work["possibility_flag"].eq(1)
            & work["exists_flag"].eq(0)
        ].copy()
        g_zero = _filter_window_fit_absent(g_zero, g_active, case, spec["pair"])
        g_active = g_active.replace([np.inf, -np.inf], np.nan).dropna(subset=[xvar, yvar, "E2_present"])
        if xvar == "E2_present" or yvar == "E2_present":
            g_zero = g_zero.assign(E2_present=0.0)
        g_zero = g_zero.replace([np.inf, -np.inf], np.nan).dropna(subset=[xvar, yvar])
        return g_active, g_zero

    def _plot_geometry_bubble(source_df, scope_label, scope_slug, case, spec, out_dir=FIGPRESENTDIR):
        filename = f"FIG_GEOM_BUBBLE_{scope_slug}_{_geom_case_slug(case)}.png"

        def _save_current():
            path = os.path.join(out_dir, filename)
            plt.savefig(path, dpi=230, bbox_inches="tight", facecolor=plt.gcf().get_facecolor())
            plt.close()
            return path

        def _empty_plot(reason):
            fig, ax = plt.subplots(figsize=(9.2, 7.0), facecolor="white")
            ax.axis("off")
            ax.set_title(f"{case} geometry bubble plot — {scope_label}", fontsize=14, fontweight="bold", pad=13)
            ax.text(0.5, 0.5, reason, transform=ax.transAxes, ha="center", va="center", fontsize=12)
            return _save_current()

        if spec["pair"] is None:
            return _empty_plot("No plottable interactions with available coordinates.")
        xvar, yvar = spec["pair"]
        g_active, g_zero = _plot_bits(source_df, case, spec)
        if g_active.empty and g_zero.empty:
            return _empty_plot("No active or window-fit absent interactions after coordinate availability filtering.")

        fig, ax = plt.subplots(figsize=(9.2, 7.0), facecolor="white")
        _geom_or_hist(ax, g_active, g_zero, spec, for_panel=False)
        plot_kind = "geometry-frequency histogram" if spec.get("histogram_mode") else "geometry bubble plot"
        ax.set_title(f"{case} {plot_kind} — {scope_label}", fontsize=14, fontweight="bold", pad=13)
        if spec.get("histogram_mode"):
            ax.set_xlabel(GEOM_LABELS.get(spec["hist_geom_var"], spec["hist_geom_var"]), fontsize=12, fontweight="bold")
            ax.set_ylabel("Interaction frequency", fontsize=12, fontweight="bold")
        else:
            ax.set_xlabel(GEOM_LABELS.get(xvar, xvar), fontsize=12, fontweight="bold")
            ax.set_ylabel(GEOM_LABELS.get(yvar, yvar), fontsize=12, fontweight="bold")
        ax.grid(alpha=0.25)
        ax.legend(title="Conformer", frameon=False, loc="best")
        mode_note = "Histogram bars = number of interactions per geometry range" if spec.get("histogram_mode") else "Bubble size = E2"
        ax.text(
            0.01, 0.01,
            f"{mode_note}; black dots = window-fit absent interactions; active n = {len(g_active)}; black-dot n = {len(g_zero)}",
            transform=ax.transAxes,
            ha="left", va="bottom", fontsize=8,
            bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "alpha": 0.82, "edgecolor": "#BBBBBB"},
        )
        plt.subplots_adjust(left=0.11, right=0.97, top=0.89, bottom=0.12)
        return _save_current()

    _active_geom_core = df[df["core_interaction"].eq(True)].copy()
    bubble_sources = [
        ("core interactions only (≤27)", "core", _active_geom_core),
        ("biomotifs", "biomotifs", df_biomotifs.copy()),
    ]
    # same bins and scales for the two csvs
    CASE_PLOT_SPECS = {case: _mk_plot_spec(case, [_active_geom_core, df_biomotifs.copy()]) for case in GEOM_BUBBLE_CASES}
    _geom_count_rows = []
    for _scope_label, _scope_slug, _src in bubble_sources:
        for _case in GEOM_BUBBLE_CASES:
            _spec = CASE_PLOT_SPECS[_case]
            _case_df = _filter_case_rows(_src, _case)
            _active_before = _case_df[
                _case_df["possibility_flag"].eq(1)
                & _case_df["exists_flag"].eq(1)
                & _case_df["E2_present"].notna()
            ].copy()
            _absent_before = _case_df[
                _case_df["possibility_flag"].eq(1)
                & _case_df["exists_flag"].eq(0)
            ].copy()
            _absent_window = _filter_window_fit_absent(_absent_before, _active_before, _case, _spec["pair"]) if _spec["pair"] is not None else _absent_before.iloc[0:0].copy()
            _active_tmp, _absent_tmp = _plot_bits(_src, _case, _spec) if _spec["pair"] is not None else (_src.iloc[0:0].copy(), _src.iloc[0:0].copy())
            _geom_count_rows.append({
                "scope": _scope_label,
                "paper_geometry_case": _case,
                "active_interactions_before_coordinate_filter": len(_active_before),
                "window_fit_absent_interactions_before_coordinate_filter": len(_absent_window),
                "plotted_coordinate_pair": "" if _spec["pair"] is None else (f"{_spec['hist_geom_var']} histogram vs interaction frequency" if _spec.get("histogram_mode") else f"{_spec['pair'][0]} vs {_spec['pair'][1]}"),
                "plotted_active_interactions_after_coordinate_filter": len(_active_tmp),
                "plotted_window_fit_absent_interactions_after_coordinate_filter": len(_absent_tmp),
            })
            _plot_geometry_bubble(_src, _scope_label, _scope_slug, _case, _spec)
    save_table_to(pd.DataFrame(_geom_count_rows), SUPPTABLEDIR, "Raw_Figure_Geometry_Bubble_Dataset_Counts.csv")

    def _paper_geom(ax, source_df, case, spec):
        if spec["pair"] is None:
            ax.text(0.5, 0.5, "No plottable data", transform=ax.transAxes, ha="center", va="center", fontsize=10)
            ax.set_xticks([])
            ax.set_yticks([])
            _style_paper_axes(ax)
            return
        xvar, yvar = spec["pair"]
        g_active, g_zero = _plot_bits(source_df, case, spec)
        _geom_or_hist(ax, g_active, g_zero, spec, for_panel=True)
        if spec.get("histogram_mode"):
            ax.set_xlabel(GEOM_LABELS.get(spec["hist_geom_var"], spec["hist_geom_var"]), fontsize=9)
            ax.set_ylabel("Interaction frequency", fontsize=9)
        else:
            ax.set_xlabel(GEOM_LABELS.get(xvar, xvar), fontsize=9)
            ax.set_ylabel(GEOM_LABELS.get(yvar, yvar), fontsize=9)
        ax.grid(alpha=0.25)
        _style_paper_axes(ax)
        handles, labels = ax.get_legend_handles_labels()
        if handles:
            ax.legend(handles, labels, loc="best", frameon=True, framealpha=0.95, fontsize=8, title=_legend_case_label(case), title_fontsize=8, borderpad=0.35)

    def _six_geom_panels(source_df, filename, saver):
        fig, axes = plt.subplots(2, 3, figsize=(18.0, 11.8), facecolor="white")
        axes = np.array(axes).reshape(-1)
        for ax, case in zip(axes, GEOM_BUBBLE_CASES):
            _paper_geom(ax, source_df, case, CASE_PLOT_SPECS[case])
        plt.subplots_adjust(left=0.06, right=0.985, top=0.985, bottom=0.07, wspace=0.22, hspace=0.24)
        saver(filename)

    _six_geom_panels(_active_geom_core, "FIG_GEOM_BUBBLE_core_6panel.png", savefig_paper)
    _six_geom_panels(df_biomotifs.copy(), "FIG_GEOM_BUBBLE_biomotifs_6panel.png", savefig_supplementary)

    STACK_COLORS = {
        "n→π*": "#3B82F6", "n→σ*": "#FF7A1A", "π→σ*": "#A3A3A3", "σ→π*": "#FACC15", "σ→σ*": "#22C55E",
    }
    NO_NPI_TYPES = [t for t in EXPECTED_MAJOR_TYPES if t != "n→π*"]

    def _stack_matrix(source_df, sub=None, types=None):
        types = types or EXPECTED_MAJOR_TYPES
        g = source_df.copy() if sub is None else source_df[source_df["substituent"].eq(sub)].copy()
        mat = g.groupby(["conformer", "type_simple"], observed=False).agg(total_E2=("total_E2", "sum")).reset_index()
        mat = mat.pivot_table(index="conformer", columns="type_simple", values="total_E2", aggfunc="sum", fill_value=0)
        return mat.reindex(index=CONFORMER_ORDER, columns=types, fill_value=0)

    def _plot_3d_stacked(mat, title, filename, zlabel, types):
        if mat.empty:
            return None
        mat = mat.reindex(columns=types).fillna(0).astype(float)
        fig = plt.figure(figsize=(13.6, 8.2), facecolor="#F1F1F1")
        ax = fig.add_subplot(111, projection="3d")
        ax.set_facecolor("#F1F1F1")

        x = np.arange(len(mat.index), dtype=float)
        dx, dy = 0.56, 0.82
        y0 = -dy / 2
        bottoms = np.zeros(len(x), dtype=float)
        totals = mat.sum(axis=1).to_numpy(dtype=float)
        ymax = max(float(totals.max()), 0.1)
        thin_threshold = max(0.11 * ymax, 0.06)
        label_step = max(0.060 * ymax, 0.05)
        label_y = y0 - 0.34
        overlay_specs = []

        def _project_to_figure(x3, y3, z3):
            x2, y2, _ = proj3d.proj_transform(x3, y3, z3, ax.get_proj())
            xd, yd = ax.transData.transform((x2, y2))
            return fig.transFigure.inverted().transform((xd, yd))

        for xi in x:
            shadow = [
                (xi - dx * 0.65, y0 - 0.08, 0),
                (xi + dx * 0.60, y0 - 0.08, 0),
                (xi + dx * 1.05, y0 + dy * 0.82, 0),
                (xi - dx * 0.20, y0 + dy * 0.82, 0),
            ]
            ax.add_collection3d(Poly3DCollection([shadow], facecolors=(0, 0, 0, 0.11), edgecolors="none"))

        outside_z = np.zeros(len(x), dtype=float)
        thin_counts = np.zeros(len(x), dtype=int)
        for typ in types:
            vals = mat[typ].to_numpy(dtype=float) if typ in mat.columns else np.zeros(len(x))
            ax.bar3d(
                x - dx / 2, y0, bottoms, dx, dy, vals,
                color=STACK_COLORS.get(typ, "#DDDDDD"), edgecolor="#55606A",
                linewidth=0.42, shade=True, zsort="average"
            )
            for i, val in enumerate(vals):
                if val <= 0:
                    continue
                top_z = bottoms[i] + val
                if val >= thin_threshold:
                    overlay_specs.append({
                        "kind": "inside",
                        "bar": i,
                        "anchor": (x[i], y0 - 0.015, bottoms[i] + val / 2),
                        "label": f"{val:.2f}",
                        "fontsize": 9.6,
                    })
                else:
                    thin_counts[i] += 1
                    label_z = max(top_z + 0.028 * ymax, outside_z[i] + label_step)
                    outside_z[i] = label_z
                    overlay_specs.append({
                        "kind": "outside",
                        "bar": i,
                        "level": int(thin_counts[i]),
                        "anchor": (x[i], y0 + 0.03, top_z),
                        "elbow": (x[i], label_y + 0.11, label_z),
                        "label_base": (x[i], label_y, label_z),
                        "label": f"{val:.2f}",
                        "fontsize": 9.6,
                    })
            bottoms += vals

        for i, total in enumerate(totals):
            overlay_specs.append({
                "kind": "total",
                "bar": i,
                "anchor": (x[i], y0 + dy * 0.28, total + 0.060 * ymax),
                "label": f"Σ {total:.2f}",
                "fontsize": 10.6,
            })

        ax.set_title(title, fontsize=16, fontweight="bold", pad=18, color="#142632")
        ax.set_zlabel(zlabel, fontsize=12, fontweight="bold", labelpad=10)
        ax.set_xticks(x)
        ax.set_xticklabels(list(mat.index), fontsize=12, fontweight="bold")
        ax.set_yticks([])
        ax.set(xlim=(-0.55, len(x) - 0.45), ylim=(label_y - 0.20, y0 + dy + 0.35), zlim=(0, ymax * 1.34))
        ax.view_init(elev=16, azim=-76)
        try:
            ax.set_box_aspect((1.70, 0.82, 1.12))
        except Exception:
            pass
        for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
            axis.pane.set_facecolor((0.955, 0.955, 0.955, 1))
            axis.pane.set_edgecolor("#CACFD3")
        ax.zaxis._axinfo["grid"].update(color="#BEC3C7", linewidth=0.80, linestyle="-")
        ax.xaxis._axinfo["grid"]["linewidth"] = 0
        ax.yaxis._axinfo["grid"]["linewidth"] = 0
        ax.tick_params(axis="z", labelsize=10, pad=4)

        handles = [plt.Rectangle((0, 0), 1, 1, facecolor=STACK_COLORS.get(t, "#DDDDDD"), edgecolor="#55606A", linewidth=0.45) for t in types]
        fig.legend(handles, types, loc="lower center", bbox_to_anchor=(0.5, 0.025), ncol=max(1, len(types)), frameon=False, fontsize=10.3, handlelength=1.15, handletextpad=0.40, columnspacing=0.9)
        plt.subplots_adjust(left=0.03, right=0.97, top=0.90, bottom=0.15)

        fig.canvas.draw()
        overlay_ax = fig.add_axes([0, 0, 1, 1], frameon=False, zorder=100)
        overlay_ax.set_axis_off()
        overlay_ax.set_xlim(0, 1)
        overlay_ax.set_ylim(0, 1)

        for spec in overlay_specs:
            if spec["kind"] in {"inside", "total"}:
                xf, yf = _project_to_figure(*spec["anchor"])
                overlay_ax.text(xf, yf, spec["label"], ha="center", va="center", fontsize=spec["fontsize"], fontweight="bold", color="black")
            else:
                axf, ayf = _project_to_figure(*spec["anchor"])
                exf, eyf = _project_to_figure(*spec["elbow"])
                lxf, lyf = _project_to_figure(*spec["label_base"])
                sign = -1 if spec["bar"] < len(x) / 2 else 1
                xoff = sign * (0.012 + 0.006 * ((spec["level"] - 1) % 2))
                yoff = 0.006 * (spec["level"] - 1)
                lxf = min(max(lxf + xoff, 0.05), 0.95)
                lyf = min(max(lyf + yoff, 0.08), 0.95)
                overlay_ax.plot([axf, exf, lxf], [ayf, eyf, lyf], color="#6F7880", linewidth=0.95)
                overlay_ax.text(lxf, lyf, spec["label"], ha="center", va="center", fontsize=spec["fontsize"], fontweight="bold", color="black")

        path = os.path.join(FIGPRESENTDIR, filename)
        plt.savefig(path, dpi=300, bbox_inches="tight", facecolor=fig.get_facecolor())
        plt.close()
        return path

    for _sub in SUB_ORDER:
        _plot_3d_stacked(_stack_matrix(type_sample, _sub, EXPECTED_MAJOR_TYPES), f"Stacked interaction-energy profile — {_sub}", f"FIG_3D_STACKED_{_sub}.png", "Interaction energy (kcal/mol)", EXPECTED_MAJOR_TYPES)
    _plot_3d_stacked(_stack_matrix(type_sample, None, EXPECTED_MAJOR_TYPES), "Stacked total interaction energy — all substituents", "FIG_3D_STACKED_total_all_substituents.png", "Total interaction energy (kcal/mol)", EXPECTED_MAJOR_TYPES)

    for _sub in SUB_ORDER:
        _plot_3d_stacked(_stack_matrix(type_sample, _sub, NO_NPI_TYPES), f"Stacked non-n→π* interaction-energy profile — {_sub}", f"FIG_3D_STACKED_no_npi_{_sub}.png", "Interaction energy excluding n→π* (kcal/mol)", NO_NPI_TYPES)
    _plot_3d_stacked(_stack_matrix(type_sample, None, NO_NPI_TYPES), "Stacked total non-n→π* interaction energy — all substituents", "FIG_3D_STACKED_no_npi_total_all_substituents.png", "Total interaction energy excluding n→π* (kcal/mol)", NO_NPI_TYPES)

    for _sub in SUB_ORDER:
        _plot_3d_stacked(_stack_matrix(type_sample_core, _sub, EXPECTED_MAJOR_TYPES), f"Stacked core-interaction energy profile — {_sub}", f"FIG_3D_STACKED_core_{_sub}.png", "Core interaction energy (kcal/mol)", EXPECTED_MAJOR_TYPES)
    _plot_3d_stacked(_stack_matrix(type_sample_core, None, EXPECTED_MAJOR_TYPES), "Stacked total core-interaction energy — all substituents", "FIG_3D_STACKED_core_total_all_substituents.png", "Total core interaction energy (kcal/mol)", EXPECTED_MAJOR_TYPES)

    for _sub in SUB_ORDER:
        _plot_3d_stacked(_stack_matrix(type_sample_core, _sub, NO_NPI_TYPES), f"Stacked core non-n→π* interaction-energy profile — {_sub}", f"FIG_3D_STACKED_core_no_npi_{_sub}.png", "Core interaction energy excluding n→π* (kcal/mol)", NO_NPI_TYPES)
    _plot_3d_stacked(_stack_matrix(type_sample_core, None, NO_NPI_TYPES), "Stacked total core non-n→π* interaction energy — all substituents", "FIG_3D_STACKED_core_no_npi_total_all_substituents.png", "Total core interaction energy excluding n→π* (kcal/mol)", NO_NPI_TYPES)

if not E_core.empty:
    same_interaction_stats = pd.DataFrame({
        "interaction_id": E_core.columns,
        "interaction_code": [code_map.get(c) for c in E_core.columns],
        "type_simple": [type_map.get(c) for c in E_core.columns],
        "long_label": [label_map.get(c) for c in E_core.columns],
        "total_E2": E_core.sum(axis=0, skipna=True).values,
        "variance_E2": E_core.var(axis=0, skipna=True).values,
        "max_E2": E_core.max(axis=0, skipna=True).values,
        "n_nonzero_samples": (E_core.fillna(0) > 0).sum(axis=0).values,
    }).sort_values(["variance_E2", "total_E2"], ascending=False)

    top_profile_ids = same_interaction_stats.head(12)["interaction_id"].tolist()
    if top_profile_ids:
        prof = E_core[top_profile_ids].copy().reindex(ordered_sample_index)
        x = np.arange(len(prof.index))
        xlabels = [f"{a}_{b}" for a, b in prof.index]
        ncols = 3
        nrows = math.ceil(len(top_profile_ids) / ncols)
        fig, axes = plt.subplots(nrows, ncols, figsize=(16, max(4 * nrows, 8)), sharex=True, sharey=False)
        axes = np.array(axes).reshape(-1)
        for ax, inter in zip(axes, top_profile_ids):
            y = prof[inter].fillna(0).values
            ax.plot(x, y, marker="o")
            ax.set_title(f"{code_map.get(inter)} | {type_map.get(inter)}")
            ax.set_ylabel("E2")
            ax.grid(alpha=0.25)
            tick_idx = np.arange(0, len(xlabels), 4)
            ax.set_xticks(tick_idx)
            ax.set_xticklabels([xlabels[i] for i in tick_idx], rotation=90, fontsize=7)
        for ax in axes[len(top_profile_ids):]:
            ax.axis("off")
        for ax in axes[:len(top_profile_ids)]:
            ax.set_xlabel("Substituent_conformer")
        fig.suptitle("Top same core interactions across all samples", y=1.02, fontsize=14)
        savefig("FIG_03B_same_core_interaction_profiles_top12.png")

    type_profile = type_sample_core.groupby(["substituent", "sigma_p", "conformer", "type_simple"], observed=False).agg(total_E2=("total_E2", "sum")).reset_index()
    type_list = type_order(type_profile["type_simple"])
    if type_list:
        ncols = 2
        nrows = math.ceil(len(type_list) / ncols)
        fig, axes = plt.subplots(nrows, ncols, figsize=(14, max(4 * nrows, 6)), sharex=True)
        axes = np.array(axes).reshape(-1)
        for ax, typ in zip(axes, type_list):
            g = type_profile[type_profile["type_simple"].eq(typ)].copy()
            for conf in CONFORMER_ORDER:
                gg = g[g["conformer"].eq(conf)].sort_values("sigma_p")
                if not gg.empty:
                    ax.plot(gg["substituent"], gg["total_E2"], marker="o", label=conf)
            ax.set_title(f"{typ} (core atoms only)")
            ax.set_ylabel("ΣE2")
            ax.grid(alpha=0.25)
        for ax in axes[len(type_list):]:
            ax.axis("off")
        for ax in axes[:len(type_list)]:
            ax.set_xlabel("Substituent, EDG → EWG")
        handles, labels = axes[0].get_legend_handles_labels()
        if handles:
            axes[0].legend(handles, labels, title="Conformer", fontsize=8)
        fig.suptitle("Core interaction-type profiles", y=1.02, fontsize=14)
        savefig("FIG_03C_core_interaction_type_profiles.png")

plot = pn_total.sort_values("sigma_p")
fig, ax = plt.subplots(figsize=(10, 5))
ax.plot(plot["substituent"], plot["N_mean_npi_E2"], marker="o", label="N mean")
ax.plot(plot["substituent"], plot["P_mean_npi_E2"], marker="o", label="P mean")
ax.set_title("P vs N n→π* stabilization")
ax.set_xlabel("Substituent, EDG → EWG")
ax.set_ylabel("Mean ΣE2(n→π*)")
ax.grid(alpha=0.25)
ax.legend()
savefig("FIG_04_P_vs_N_npi.png")

fig, axes = plt.subplots(2, 1, figsize=(11, 9), sharex=True)
for ax, src, title in [
    (axes[0], summary, "Non-n→π* stabilization across EDG→EWG substituent series (all atoms)"),
    (axes[1], summary_core, f"Non-n→π* stabilization across EDG→EWG substituent series (core atoms only, atom index ≤ {CORE_MAX_ATOM})"),
]:
    for conf in CONFORMER_ORDER:
        g = src[src["conformer"].eq(conf)].sort_values("sigma_p")
        ax.plot(g["substituent"], g["non_npi_total_E2"], marker="o", label=conf)
    ax.set_title(title)
    ax.set_ylabel("ΣE2 excluding n→π*")
    ax.legend(title="Conformer")
    ax.grid(alpha=0.25)
axes[-1].set_xlabel("Substituent, EDG → EWG")
savefig("FIG_04B_non_npi_E2_EDG_EWG_by_conformer.png")

fig, axes = plt.subplots(1, 2, figsize=(14, 5))
for conf in CONFORMER_ORDER:
    g = summary[summary["conformer"].eq(conf)].sort_values("sigma_p")
    axes[0].plot(g["sigma_p"], g["npi_total_E2"], marker="o", label=conf)
    axes[1].plot(g["sigma_p"], g["total_E2"], marker="o", label=conf)
axes[0].set_title("n→π* along EDG→EWG axis"); axes[0].set_xlabel("Hammett σp"); axes[0].set_ylabel("ΣE2(n→π*)"); axes[0].legend(); axes[0].grid(alpha=0.25)
axes[1].set_title("Total stabilization along EDG→EWG axis"); axes[1].set_xlabel("Hammett σp"); axes[1].set_ylabel("ΣE2(all)"); axes[1].legend(); axes[1].grid(alpha=0.25)
savefig("FIG_05_hammett_trends.png")

geom_top = geom_corr[geom_corr["n_points"] >= MIN_CORR_N].head(6)
if not geom_top.empty:
    nplots = len(geom_top); ncols = 3; nrows = math.ceil(nplots/ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(5*ncols, 4*nrows))
    axes = np.array(axes).reshape(-1)
    for ax, (_, r) in zip(axes, geom_top.iterrows()):
        typ, fam, var = r["type_simple"], r["conformer_family"], r["geometry_variable"]
        g = existing[(existing["type_simple"].eq(typ)) & ((existing["conformer_family"].eq(fam)) if fam != "all" else True)].copy() if fam != "all" else existing[existing["type_simple"].eq(typ)].copy()
        for cfam in ["N", "P"]:
            gg = g[g["conformer_family"].eq(cfam)]
            if not gg.empty:
                ax.scatter(gg[var], gg["E2_present"], alpha=0.75, label=cfam)
        ax.set_title(f"{typ}: E2 vs {var}\nr={r['corr_E2_vs_geometry']:.2f}, set={fam}")
        ax.set_xlabel(var); ax.set_ylabel("E2")
        ax.legend(fontsize=8)
    for ax in axes[nplots:]: ax.axis("off")
    savefig("FIG_06_geometry_E2_scatter_key_types.png")

pp = pareto.head(30).copy()
fig, axes = plt.subplots(1, 2, figsize=(15, 5))

top_bar = pp.head(15).sort_values("total_E2", ascending=True)
axes[0].barh(top_bar["interaction_code"], top_bar["total_E2"])
axes[0].set_title("Top individual interaction contributions")
axes[0].set_xlabel("Total E2")
axes[0].set_ylabel("Interaction code")

axes[1].plot(pp["rank"], pp["cumulative_fraction"], marker="o")
for frac in [0.50, 0.80, 0.90]:
    axes[1].axhline(frac, linestyle="--", linewidth=0.8, alpha=0.6)
    axes[1].text(pp["rank"].max(), frac, f" {int(frac*100)}%", va="center", fontsize=8)
axes[1].set_ylim(0, 1.05)
axes[1].set_title("Cumulative contribution curve")
axes[1].set_xlabel("Ranked interactions included")
axes[1].set_ylabel("Cumulative fraction of total E2")
axes[1].grid(alpha=0.25)

savefig("FIG_12_pareto_interaction_contributions.png")


geom3d = existing[["type_simple", "distance", "angle", "E2_present", "dihedral"]].dropna(subset=["distance", "angle", "E2_present"]).copy()
geom3d = geom3d[geom3d["type_simple"].isin(type_order(geom3d["type_simple"])[:4])].copy()
if not geom3d.empty:
    fig = plt.figure(figsize=(11, 8))
    ax = fig.add_subplot(111, projection='3d')
    use_dihedral = "dihedral" in geom3d.columns and geom3d["dihedral"].notna().any()
    color_vals = geom3d["dihedral"].fillna(0) if use_dihedral else geom3d["E2_present"]
    cbar_label = "Dihedral" if use_dihedral else "E2"
    cmap_name = "coolwarm" if use_dihedral else "viridis"

    finite_vals = pd.Series(color_vals).replace([np.inf, -np.inf], np.nan).dropna()
    if not finite_vals.empty:
        vmin = float(finite_vals.quantile(0.05))
        vmax = float(finite_vals.quantile(0.95))
        if use_dihedral and vmin < 0 < vmax:
            vmax_abs = max(abs(vmin), abs(vmax))
            norm = mcolors.TwoSlopeNorm(vmin=-vmax_abs, vcenter=0.0, vmax=vmax_abs)
        else:
            if vmax <= vmin:
                vmax = vmin + 1e-9
            norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    else:
        norm = None

    pts = ax.scatter(
        geom3d["distance"],
        geom3d["angle"].fillna(0),
        geom3d["E2_present"],
        c=color_vals,
        cmap=cmap_name,
        norm=norm,
        s=42,
        alpha=0.78,
        depthshade=True,
        linewidths=0.15,
        edgecolors="k"
    )
    ax.view_init(elev=24, azim=-58)
    ax.set_title("3D geometry-energy cloud (key interaction types)")
    ax.set_xlabel("Distance")
    ax.set_ylabel("Angle")
    ax.set_zlabel("E2")
    fig.colorbar(pts, ax=ax, shrink=0.7, pad=0.08, label=cbar_label)
    savefig("FIG_22_3D_geometry_energy_cloud.png")

FIGURE_TABLE_INDEX = []

def save_figure_table(df_in, filename, figure_name, purpose):
    out = df_in.copy() if isinstance(df_in, pd.DataFrame) else pd.DataFrame()
    if out.empty:
        out = pd.DataFrame({"note": ["No data available for this figure in the current run."]})
    out = tidy_round(out)
    raw_name = "Raw_Figure_" + filename
    save_table_to(out, SUPPTABLEDIR, raw_name)
    FIGURE_TABLE_INDEX.append({
        "raw_supplementary_table": raw_name,
        "figure_file": figure_name,
        "purpose": purpose,
        "n_rows": len(out),
        "n_columns": len(out.columns),
    })

def _select_cols(df_in, cols):
    return df_in[[c for c in cols if c in df_in.columns]].copy() if isinstance(df_in, pd.DataFrame) else pd.DataFrame()

save_figure_table(
    _select_cols(fig01_all.sort_values(["sigma_p", "conformer_family", "type_simple"]),
                 ["substituent", "sigma_p", "conformer_family", "type_simple", "mean_total_E2"]),
    "Table_01_FIG_01_comparable_energies_N_vs_P.csv",
    "FIG_01_reference_style_comparable_energies_N_vs_P.png",
    "Mean total E2 by substituent, conformer family, and interaction type."
)

fig01b_cols = ["substituent", "sigma_p", "conformer", "conformer_family", "type_simple", "total_E2", "count_existing", "mean_E2", "max_E2"]
save_figure_table(_select_cols(type_sample.sort_values(["sigma_p", "conformer", "type_simple"]), fig01b_cols),
                  "Table_02_FIG_01B_E2_by_substituent_type_lines_each_conformer.csv",
                  "FIG_01B_E2_by_substituent_type_lines_each_conformer.png",
                  "All interaction types; all atoms; source table for the main FIG_01B panel plot.")
save_figure_table(_select_cols(type_sample[~type_sample["type_simple"].isin(["n→σ*", "σ→σ*"])].sort_values(["sigma_p", "conformer", "type_simple"]), fig01b_cols),
                  "Table_03_FIG_01B2_no_nsigma_ssigma.csv",
                  "FIG_01B2_E2_by_substituent_type_lines_each_conformer_no_nsig_ssig.png",
                  "FIG_01B variant excluding n→σ* and σ→σ*.")
save_figure_table(_select_cols(type_sample[type_sample["type_simple"].ne("n→π*")].sort_values(["sigma_p", "conformer", "type_simple"]), fig01b_cols),
                  "Table_04_FIG_01B3_no_npi.csv",
                  "FIG_01B3_E2_by_substituent_type_lines_each_conformer_no_npi.png",
                  "FIG_01B variant excluding n→π*.")
save_figure_table(_select_cols(type_sample_core.sort_values(["sigma_p", "conformer", "type_simple"]), fig01b_cols),
                  "Table_05_FIG_01B4_core_atoms.csv",
                  "FIG_01B4_E2_by_substituent_type_lines_each_conformer_core_atoms.png",
                  f"Core-atoms-only FIG_01B source table, atom index ≤ {CORE_MAX_ATOM}.")
save_figure_table(_select_cols(type_sample_core[~type_sample_core["type_simple"].isin(["n→σ*", "σ→σ*"])].sort_values(["sigma_p", "conformer", "type_simple"]), fig01b_cols),
                  "Table_06_FIG_01B5_core_no_nsigma_ssigma.csv",
                  "FIG_01B5_E2_by_substituent_type_lines_each_conformer_core_atoms_no_nsig_ssig.png",
                  f"Core-atoms-only FIG_01B variant excluding n→σ* and σ→σ*, atom index ≤ {CORE_MAX_ATOM}.")
save_figure_table(_select_cols(type_sample_core[type_sample_core["type_simple"].ne("n→π*")].sort_values(["sigma_p", "conformer", "type_simple"]), fig01b_cols),
                  "Table_07_FIG_01B6_core_no_npi.csv",
                  "FIG_01B6_E2_by_substituent_type_lines_each_conformer_core_atoms_no_npi.png",
                  f"Core-atoms-only FIG_01B variant excluding n→π*, atom index ≤ {CORE_MAX_ATOM}.")

save_figure_table(_select_cols(type_conf.sort_values(["conformer", "type_simple"]),
                                ["conformer", "conformer_family", "type_simple", "total_E2", "count_existing", "mean_E2", "max_E2", "fraction_of_conformer_E2"]),
                  "Table_08_FIG_02_interaction_type_by_conformer.csv",
                  "FIG_02_interaction_type_by_conformer.png",
                  "Total interaction-type stabilization by conformer, all atoms.")
type_conf_core_for_table = type_sample_core.groupby(["conformer", "conformer_family", "type_simple"], observed=False).agg(
    total_E2=("total_E2", "sum"),
    count_existing=("count_existing", "sum"),
    mean_E2=("mean_E2", "mean"),
    max_E2=("max_E2", "max"),
).reset_index()
type_conf_core_for_table["fraction_of_conformer_E2"] = type_conf_core_for_table["total_E2"] / type_conf_core_for_table.groupby("conformer")["total_E2"].transform("sum").replace(0, np.nan)
save_figure_table(_select_cols(type_conf_core_for_table.sort_values(["conformer", "type_simple"]),
                                ["conformer", "conformer_family", "type_simple", "total_E2", "count_existing", "mean_E2", "max_E2", "fraction_of_conformer_E2"]),
                  "Table_09_FIG_02A_interaction_type_by_conformer_core_atoms.csv",
                  "FIG_02A_interaction_type_by_conformer_core_atoms.png",
                  f"Total interaction-type stabilization by conformer for core atoms only, atom index ≤ {CORE_MAX_ATOM}.")

save_figure_table(_select_cols(globals().get("same_interaction_stats", pd.DataFrame()).head(12),
                                ["interaction_code", "type_simple", "long_label", "total_E2", "variance_E2", "max_E2", "n_nonzero_samples"]),
                  "Table_11_FIG_03B_same_core_interaction_profiles_top12.csv",
                  "FIG_03B_same_core_interaction_profiles_top12.png",
                  "Top same-core interactions used in the FIG_03B profile panels.")

def _label_groups(label):
    groups = []
    for chunk in re.split(r"[;,]+", sstr(label)):
        atoms = atom_tokens(chunk)
        if atoms:
            groups.append("/".join(atoms))
        elif chunk.strip():
            groups.append(chunk.strip())
    return groups

def _unique_atom_group(series):
    vals, seen = [], set()
    for item in pd.Series(series).dropna().astype(str):
        for group in _label_groups(item):
            if group not in seen:
                seen.add(group)
                vals.append(group)
    return ", ".join(vals)

def _fig03c_core_source_table():
    base = globals().get("type_profile", type_sample_core).sort_values(["sigma_p", "conformer", "type_simple"]).copy()
    active_core = df[
        df["core_interaction"].eq(True)
        & df["possibility_flag"].eq(1)
        & df["exists_flag"].eq(1)
    ].copy()
    if active_core.empty:
        base["donor_groups"] = ""
        base["acceptor_groups"] = ""
        return base
    groups = (
        active_core.groupby(["substituent", "sigma_p", "conformer", "type_simple"], observed=False)
        .agg(donor_groups=("donor_labels", _unique_atom_group),
             acceptor_groups=("acceptor_labels", _unique_atom_group))
        .reset_index()
    )
    out = base.merge(groups, on=["substituent", "sigma_p", "conformer", "type_simple"], how="left")
    out[["donor_groups", "acceptor_groups"]] = out[["donor_groups", "acceptor_groups"]].fillna("")
    return out

save_figure_table(_select_cols(_fig03c_core_source_table(),
                                ["substituent", "sigma_p", "conformer", "type_simple", "total_E2", "donor_groups", "acceptor_groups"]),
                  "Table_12_FIG_03C_core_interaction_type_profiles.csv",
                  "FIG_03C_core_interaction_type_profiles.png",
                  "Core interaction-type profile source data with unique donor and acceptor atom groups.")

save_figure_table(_select_cols(pn_total.sort_values("sigma_p"),
                                ["substituent", "sigma_p", "N_mean_npi_E2", "P_mean_npi_E2", "P_minus_N_npi_E2"]),
                  "Table_13_FIG_04_P_vs_N_npi.csv",
                  "FIG_04_P_vs_N_npi.png",
                  "P vs N n→π* stabilization source data.")
fig04b_all = summary[["substituent", "sigma_p", "conformer", "non_npi_total_E2"]].copy(); fig04b_all["atom_scope"] = "all atoms"
fig04b_core = summary_core[["substituent", "sigma_p", "conformer", "non_npi_total_E2"]].copy(); fig04b_core["atom_scope"] = f"core atoms ≤ {CORE_MAX_ATOM}"
save_figure_table(tidy_round(pd.concat([fig04b_all, fig04b_core], ignore_index=True).sort_values(["atom_scope", "sigma_p", "conformer"])),
                  "Table_14_FIG_04B_non_npi_E2_EDG_EWG_by_conformer.csv",
                  "FIG_04B_non_npi_E2_EDG_EWG_by_conformer.png",
                  "Non-n→π* E2 by substituent and conformer, all atoms and core-only.")
save_figure_table(_select_cols(summary.sort_values(["sigma_p", "conformer"]),
                                ["substituent", "sigma_p", "conformer", "npi_total_E2", "total_E2", "non_npi_total_E2"]),
                  "Table_15_FIG_05_hammett_trends.csv",
                  "FIG_05_hammett_trends.png",
                  "n→π*, non-n→π*, and total E2 across the EDG→EWG Hammett axis.")

geom_top_for_table = geom_corr[geom_corr["n_points"] >= MIN_CORR_N].head(6).copy()
fig06_rows = []
for _, _r in geom_top_for_table.iterrows():
    _typ, _fam, _var = _r["type_simple"], _r["conformer_family"], _r["geometry_variable"]
    _g = existing[existing["type_simple"].eq(_typ)].copy()
    if _fam != "all":
        _g = _g[_g["conformer_family"].eq(_fam)].copy()
    for _, _row in _g.iterrows():
        fig06_rows.append({
            "panel_type": _typ,
            "panel_family": _fam,
            "geometry_variable": _var,
            "substituent": _row.get("substituent"),
            "conformer": _row.get("conformer"),
            "interaction_code": _row.get("interaction_code"),
            "geometry_value": _row.get(_var),
            "E2_present": _row.get("E2_present"),
        })
save_figure_table(pd.DataFrame(fig06_rows),
                  "Table_16_FIG_06_geometry_E2_scatter_key_types.csv",
                  "FIG_06_geometry_E2_scatter_key_types.png",
                  "Raw scatter points behind the plotted top geometry–E2 panels.")

save_figure_table(_select_cols(globals().get("pp", pareto.head(30)).copy(),
                                ["rank", "interaction_code", "type_simple", "long_label", "total_E2", "cumulative_E2", "cumulative_fraction"]),
                  "Table_25_FIG_12_pareto_interaction_contributions.csv",
                  "FIG_12_pareto_interaction_contributions.png",
                  "Top individual interaction contributions and cumulative fraction used in Pareto figure.")
save_figure_table(_select_cols(geom3d,
                                ["type_simple", "distance", "angle", "E2_present", "dihedral"]),
                  "Table_27_FIG_22_3D_geometry_energy_cloud.csv",
                  "FIG_22_3D_geometry_energy_cloud.png",
                  "Raw geometry-energy cloud source data.")

save_table_to(pd.DataFrame(FIGURE_TABLE_INDEX), SUPPTABLEDIR, "Raw_Figure_Table_Index.csv")

export_excel_playground(df, os.path.join(OUTDIR, EXCEL_PLAYGROUND_NAME))


print("DONE")
print(f"Output folder: {OUTDIR}")

