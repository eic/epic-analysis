import os
import os
import csv
import sys
from pathlib import Path
import ROOT
import uproot

# Function to create necessary directory structure
def create_directory_structure(base_path, path_parts):
    dir_path = Path(base_path, *path_parts)
    dir_path.mkdir(parents=True, exist_ok=True)
    return dir_path

# Function to load CSV data into a list
def load_csv_data(csv_file):
    data = []
    if csv_file.exists():
        with open(csv_file, mode='r') as csvf:
            reader = csv.reader(csvf)
            data = [row for row in reader]
    return data

# Function to check if a file exists in the loaded CSV data
def file_exists_in_csv(file_name, csv_data):
    return any(file_name in row for row in csv_data)

# Function to append data to CSV
def append_to_csv(csv_file, file_name, num_events):
    with open(csv_file, mode='a', newline='') as csvf:
        writer = csv.writer(csvf)
        writer.writerow([file_name, num_events])

# Function to count events in ROOT files
def count_events(file_path):
    ROOT.gErrorIgnoreLevel = ROOT.kFatal
    try:
        tfile = ROOT.TFile.Open(file_path)
        ttree = tfile.Get("events")
        nevents = ttree.GetEntries()
        return nevents
    except:
        return 0

# Function to display a simple progress bar
def print_progress_bar(iteration, total, prefix='', suffix='', length=50, fill='='):
    percent = ("{0:.1f}").format(100 * (iteration / float(total)))
    filled_length = int(length * iteration // total)
    bar = fill * filled_length + '-' * (length - filled_length)
    print(f'{prefix} |{bar}| {percent}% {suffix}', end='\r')
    # Print New Line on Complete
    if iteration == total: 
        print()

# Main function
def main(campaign, detector, energy):
    s3_prefix = "s3https://eics3.sdcc.bnl.gov:9000/eictest"
    xrootd_prefix = "root://dtn-eic.jlab.org//work/eic2/"
    base_path = "hpc/nevents_databases"

    dirname = f"/EPIC/RECO/{campaign}/{detector}/DIS/NC/{energy}/"
    minQ2s = [file.GetName() for file in ROOT.TSystemDirectory(xrootd_prefix + dirname, xrootd_prefix + dirname).GetListOfFiles()]

    for minQ2 in minQ2s:
        full_dirname = dirname + minQ2
        d = ROOT.TSystemDirectory(xrootd_prefix + full_dirname, xrootd_prefix + full_dirname)

        processed_files = 0
        total_files = len(d.GetListOfFiles())
        print("Found",total_files,"for",minQ2)
        for ifile,file in enumerate(d.GetListOfFiles()):
            file_path = xrootd_prefix + full_dirname + "/" + file.GetName()
            col1_value = s3_prefix + full_dirname + "/" + file.GetName()
            path_parts = col1_value.split('/')[5:12]  # Adjust indices as needed
            dir_path = create_directory_structure(base_path, path_parts)
            csv_file = dir_path / 'data.csv'
            if ifile==0:
                csv_data = load_csv_data(csv_file)
            
            processed_files += 1
            print_progress_bar(processed_files, total_files, prefix='Processing [{}]:'.format(minQ2), suffix='Complete', length=50)
            if file_exists_in_csv(col1_value, csv_data):
                continue
            else:
                nevents = count_events(file_path)
                col2_value = nevents
                append_to_csv(csv_file, col1_value, col2_value)
        
    
if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python3 count_events.py <campaign> <detector> <energy>")
        sys.exit(1)

    campaign = sys.argv[1]
    if "epic." in campaign:
        campaign=campaign[5:] # Refactor 
    detector = sys.argv[2]
    energy = sys.argv[3]
    print("Running count_events.py for CAMPAIGN: {} ... ENERGY: {}".format(campaign,energy))
    main(campaign, detector, energy)
