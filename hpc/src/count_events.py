import os
import csv
import sys
import subprocess
from pathlib import Path
import ROOT

# --------------------------- CONFIGURATION --------------------------- #
XROOTD_SERVER = "root://dtn-eic.jlab.org"
VOLATILE_PREFIX = "/volatile/eic"
BASE_PATH = Path("hpc/nevents_databases")
S3_PREFIX = "s3https://eics3.sdcc.bnl.gov:9000/eictest"
# --------------------------------------------------------------------- #


def create_directory_structure(base_path, path_parts):
    dir_path = Path(base_path, *path_parts)
    dir_path.mkdir(parents=True, exist_ok=True)
    return dir_path


def load_csv_data(csv_file):
    if csv_file.exists():
        with open(csv_file, mode='r') as f:
            return [row for row in csv.reader(f)]
    return []


def file_exists_in_csv(file_name, csv_data):
    return any(file_name == row[0] for row in csv_data)


def append_to_csv(csv_file, file_name, num_events):
    with open(csv_file, mode='a', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([file_name, num_events])


def count_events(file_path):
    """Count entries in the 'events' TTree."""
    ROOT.gErrorIgnoreLevel = ROOT.kFatal
    try:
        tfile = ROOT.TFile.Open(file_path)
        if not tfile or tfile.IsZombie():
            return 0
        ttree = tfile.Get("events")
        if not ttree:
            return 0
        nevents = int(ttree.GetEntries())
        tfile.Close()
        return nevents
    except Exception:
        return 0


def print_progress_bar(iteration, total, prefix='', suffix='', length=50, fill='='):
    percent = "{0:.1f}".format(100 * (iteration / float(total)))
    filled_len = int(length * iteration // total)
    bar = fill * filled_len + '-' * (length - filled_len)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end='', flush=True)
    if iteration == total:
        print()


def list_remote_dirs(server, path):
    """Return subdirectories under an XRootD directory."""
    try:
        result = subprocess.check_output(
            ["xrdfs", server, "ls", path], text=True
        ).strip().splitlines()
        # keep only directory paths (heuristic)
        return [p for p in result if "/q2_" in p or "/minQ2" in p]
    except subprocess.CalledProcessError:
        return []


def list_remote_files(server, path):
    """Return .root files under an XRootD directory."""
    try:
        result = subprocess.check_output(
            ["xrdfs", server, "ls", path], text=True
        ).strip().splitlines()
        return [p for p in result if p.endswith(".root")]
    except subprocess.CalledProcessError:
        return []


def main(campaign, detector, energy):
    # Determine base path for this dataset
    if "10x166" in energy:
        data_dir = f"{VOLATILE_PREFIX}/EPIC/RECO/{campaign}/{detector}/DIS/BeAGLE1.03.02-1.0/eHe3/{energy}"
    else:
        data_dir = f"{VOLATILE_PREFIX}/EPIC/RECO/{campaign}/{detector}/DIS/NC/{energy}"

    print(f"Scanning {data_dir} ...")

    q2_dirs = list_remote_dirs(XROOTD_SERVER, data_dir)
    if not q2_dirs:
        print(f"No Q² directories found under {data_dir}")
        sys.exit(1)

    for q2_dir in q2_dirs:
        q2_label = Path(q2_dir).name
        root_files = list_remote_files(XROOTD_SERVER, q2_dir)
        total_files = len(root_files)
        if total_files == 0:
            print(f"No ROOT files found in {q2_label}")
            continue

        print(f"Found {total_files} ROOT files for {q2_label}")
        csv_data_cache = {}

        for i, file_path in enumerate(root_files, 1):
            basename = Path(file_path).name
            col1_value = S3_PREFIX + file_path.replace(VOLATILE_PREFIX, "")
            path_parts = file_path.split("/")[5:12]  # e.g. EPIC/.../q2_100to1000
            dir_path = create_directory_structure(BASE_PATH, path_parts)
            csv_file = dir_path / "data.csv"

            if csv_file not in csv_data_cache:
                csv_data_cache[csv_file] = load_csv_data(csv_file)

            csv_data = csv_data_cache[csv_file]
            print_progress_bar(i, total_files, prefix=f'[{q2_label}]', suffix='Complete')

            if file_exists_in_csv(col1_value, csv_data):
                continue

            full_remote_path = f"{XROOTD_SERVER}/{file_path}"
            nevents = count_events(full_remote_path)
            append_to_csv(csv_file, col1_value, nevents)
            csv_data.append([col1_value, nevents])

        print(f"Finished {q2_label}")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python3 count_events.py <campaign> <detector> <energy>")
        sys.exit(1)

    campaign = sys.argv[1].replace("epic.", "")  # strip 'epic.' prefix
    detector = sys.argv[2]
    energy = sys.argv[3]
    print(f"Running count_events.py for CAMPAIGN: {campaign} ... ENERGY: {energy}")
    main(campaign, detector, energy)
