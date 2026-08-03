from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
INPUT_DIR = ROOT / "input"
PDB_DIR = ROOT / "pdb"
OUTPUT_DIR = ROOT / "output"
OUTPUT_DIR.mkdir(exist_ok=True)

TOOLS = {}


def command(name):
    def decorator(func):
        TOOLS[name] = func
        return func
    return decorator


def run_script(script_name, *args):
    script = ROOT / script_name
    if not script.exists():
        print(f"ERROR: {script_name} not found")
        return False
    result = subprocess.run(
        [sys.executable, str(script), *map(str, args)],
        cwd=ROOT,
    )
    return result.returncode == 0


def move_to_output(path):
    source = Path(path)
    if not source.exists():
        return None
    destination = OUTPUT_DIR / source.name
    if destination.exists():
        destination.unlink()
    shutil.move(str(source), str(destination))
    return destination


@command("help")
def help_command():
    print("\nAvailable commands:\n")
    for name in sorted(TOOLS):
        print(f"  {name}")
    print()


@command("xyz")
def xyz():
    logs = sorted(INPUT_DIR.glob("*.log"))
    if not logs:
        print("No log files found")
        return
    run_script("molecule_visualizer.py", *logs)


@command("parse")
def parse():
    logs = sorted(INPUT_DIR.glob("*.log"))
    if not logs:
        print("No log files found")
        return
    run_script("parser.py", *logs)


@command("filter")
def filter_cmd():
    csv_files = sorted(INPUT_DIR.glob("*.csv"))
    if not csv_files:
        print("No CSV files found")
        return
    for csv_file in csv_files:
        xyz_file = INPUT_DIR / f"{csv_file.stem}.xyz"
        pdb_file = PDB_DIR / f"{csv_file.stem}.pdb"
        structure = xyz_file if xyz_file.exists() else pdb_file if pdb_file.exists() else None
        if structure is None:
            print(f"No structure found for {csv_file.name}")
            continue
        if run_script("filter.py", csv_file, structure):
            move_to_output(ROOT / f"filtered_{csv_file.stem}.csv")


@command("classify")
def classify():
    files = sorted(OUTPUT_DIR.glob("filtered_*.csv"))
    if not files:
        print("No filtered files found")
        return
    run_script("classifier.py", *files)


def run_pipeline(use_pdb=False):
    logs = sorted(INPUT_DIR.glob("*.log"))
    if not logs:
        print("No logs found")
        return

    print(f"\nFound {len(logs)} log files")

    print("\n[1/4] Creating XYZ")
    if not run_script("molecule_visualizer.py", *logs):
        print("Pipeline stopped during XYZ creation")
        return

    print("\n[2/4] Parsing")
    if not run_script("parser.py", *logs):
        print("Pipeline stopped during parsing")
        return

    print("\n[3/4] Filtering")
    filtered_files = []
    for log in logs:
        csv_file = INPUT_DIR / f"{log.stem}.csv"
        xyz_file = INPUT_DIR / f"{log.stem}.xyz"
        pdb_file = PDB_DIR / f"{log.stem}.pdb"
        structure = pdb_file if use_pdb else xyz_file

        if not csv_file.exists() or not structure.exists():
            print(f"Missing input for {log.stem}")
            return
        if not run_script("filter.py", csv_file, structure):
            print(f"Pipeline stopped while filtering {log.stem}")
            return

        filtered = move_to_output(ROOT / f"filtered_{log.stem}.csv")
        if filtered is None:
            print(f"Filtered output missing for {log.stem}")
            return
        filtered_files.append(filtered)

    print("\n[4/4] Classifying")
    if not run_script("classifier.py", *filtered_files):
        print("Pipeline stopped during classification")
        return

    print(f"\nPipeline complete -> output/ ({len(logs)} files)")


@command("pipeline")
def pipeline():
    run_pipeline(use_pdb=False)


@command("pipeline-pdb")
def pipeline_pdb():
    run_pipeline(use_pdb=True)


@command("quit")
def quit_command():
    raise SystemExit


def main():
    print("\nChemTools")
    print("Type 'help' for commands.\n")
    while True:
        try:
            name = input("chemtools> ").strip()
            if not name:
                continue
            if name not in TOOLS:
                print(f"Unknown command '{name}'")
                continue
            TOOLS[name]()
        except KeyboardInterrupt:
            print()
        except SystemExit:
            break
        except Exception as exc:
            print(f"ERROR: {exc}")


if __name__ == "__main__":
    main()
