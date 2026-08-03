import sys
from pathlib import Path

ELEMENTS = [
    "", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
    "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
]


def extract_coords(text):
    lines = text.splitlines()
    natoms = next((int(line.split()[1]) for line in lines if "NAtoms=" in line), None)
    if natoms is None:
        return None
    #The last Gaussian orientation before the main job terminates is the geometry used by the corresponding final NBO analysis
    end = next(
        (i for i, line in enumerate(lines) if "Normal termination of Gaussian" in line),
        len(lines),
    )
    headers = [
        i for i, line in enumerate(lines[:end])
        if "Input orientation:" in line or "Standard orientation:" in line
    ]
    if not headers:
        return None
    coords = []
    for line in lines[headers[-1] + 5: headers[-1] + 5 + natoms]:
        parts = line.split()
        if len(parts) < 6:
            return None
        atomic_number = int(parts[1])
        if not 0 < atomic_number < len(ELEMENTS):
            return None
        coords.append((ELEMENTS[atomic_number], *map(float, parts[3:6])))
    return coords if len(coords) == natoms else None


def write_xyz(coords, out_path):
    path = Path(out_path).with_suffix(".xyz")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write(f"{len(coords)}\nGaussian geometry\n")
        for atom, x, y, z in coords:
            handle.write(f"{atom} {x:.8f} {y:.8f} {z:.8f}\n")
    return str(path)


def main():
    if len(sys.argv) < 2:
        print("Usage: python molecule_visualizer.py file1.log [file2.log ...]")
        return 1
    failed = False
    for filename in sys.argv[1:]:
        path = Path(filename)
        if not path.exists():
            print("File not found:", path)
            failed = True
            continue
        coords = extract_coords(path.read_text(encoding="utf-8", errors="ignore"))
        if not coords:
            print(f"No usable Gaussian geometry found in {path.name}")
            failed = True
            continue
        print(write_xyz(coords, path.with_suffix(".xyz")))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
