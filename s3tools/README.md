# S3 Tools

Scripts to provide access to S3 and interface with this SIDIS analysis
repository

For more details, see [S3 file storage documentation](https://doc.athena-eic.org/en/latest/howto/s3_file_storage.html)

## Quick Start

### Full Simulations
One top-level script automates all the work:
- `s3tools/make-canyonlands-config.sh` (best to run from top-level directory)
  - running with no arguments will print the usage guide
  - output:
    - config file, with file names, Q2 minima, and cross sections; used as
      input to the analysis macros
    - the downloaded full simulation files, if you chose to download from S3
  - this script contains some settings such as directory paths to data on S3,
    for the convenience of SIDIS full simulation analysis
  - you may have to run `add-host.sh` first, if you have not yet; set env vars
    `$S3_ACCESS_KEY` and `$S3_SECRET_KEY` with the login and password beforehand:
    - `export=S3_ACCESS_KEY=*****`
    - `export=S3_SECRET_KEY=*****`
  - this script calls several other scripts in this directory; read on for their
    documentation
- if you are not running in the Singularity/Docker container, download and install
  the [MinIO client](https://docs.min.io/docs/minio-client-complete-guide) first
  - (otherwise update your Singularity image)

### Fast Simulations
Two options:
- download `hepmc` files from S3 and run Delphes locally
  - this is automated by `s3tools/make-fastsim-S3-config.sh`; run with no
    arguments to print the usage guide, and note the arguments are a bit
    different than those for the full simulation scripts
  - this script can download `hepmc` files from S3, run Delphes on them, and
    generate the config file
    - Delphes will run using one thread per Q2 minimum; edit the script
      to change this behavior
  - similar to full simulations, you need to have the S3 access and secret
    environment variables set in order to download `hepmc` files
- alternatively, if you already have a directory of Delphes output ROOT files,
  use use `make-fastsim-config.sh` to create a config file

## Accessing S3 Files
- first, download the [MinIO client](https://docs.min.io/docs/minio-client-complete-guide)
  - if you are using the Singularity or Docker container, it is already installed
  - the main command is `mc`
- next, setup your client to connect to our S3 host (ask someone for credentials):
  - first set env vars `$S3_ACCESS_KEY` and `$S3_SECRET_KEY` with the login and password:
    - `export=S3_ACCESS_KEY=*****`
    - `export=S3_SECRET_KEY=*****`
  - then add our S3 host to MinIO client: `add-host.sh`
    - this only needs to be done once on your machine or container
- find a directory on S3 with full simulation files; example S3 navigation commands:
  - top-level ATHENA directory list: `mc ls S3/eictest/ATHENA`
  - show directory tree: `mc tree S3/eictest/ATHENA/RECO/acadia-v2.1/DIS/NC/`
  - list files: `mc ls S3/eictest/ATHENA/RECO/acadia-v2.1/DIS/NC/10x275/minQ2=1`
- from here, you have two options, both of which are described in the next sections:
  - generate a list of S3 file URLs, which will be "streamed" when running
    the analysis code, and will not be stored locally
  - if streaming files from S3 is unsuitable, you can download them instead
    using MinIO client: `mc cp S3/.../.../source.root ./your/data/directory/`;
    see below for a download script for automation and filtering

## Generating Config Files
Next we need to make a "config file", which consists of the file name, and
additional columns such as cross section and Q2min. Follow the next sections,
whether you plan to stream from S3 or download.

### Config File Format
The config files require the following columns, in this order:
- file name (relative to the top-level directory, unless you use an absolute
  path)
- the number of events for the weighting and cross section
  - set to `0` for all
  - this is not related to `Analysis::maxEvents`, which limits how
    many events to process
- cross section (can be obtained from Pythia output logs, for example)
- minimum Q2

**Patch**: the above format is the original format, however, the current Q2min
weighting implementation requires a new format. In case we revert to using the
above old format, we temporarily use the script `reformat-config.sh` to
transform the above old format into the new format. See comments in
`reformat-config.sh` for details. Execute:
  - `s3tools/reformat-config.sh files.config files.new.config`

### Cross Sections
- cross sections are stored in `datarec/xsec/xsec.dat`; use `read-xsec-table.sh`
  to get the cross section for a particular beam energy setting and Q2 minimum
- in case you want to update `xsec.dat`:
  - the script `get-cross-section.sh` will read the cross section from
    `GenCrossSection` in a `hepmc` file; use `get-cross-section-ALL.sh` to
    automate running `get-cross-section.sh` over all `hepmc` files in a specific
    directory tree; this will populate `datarec/xsec/*.xsec` files, one for each
    `hepmc` file
  - next use `tabulate-cross-section.py` to read the tree of `.xsec` files into
    a table of cross sections, output to `datarec/xsec/xsec.dat`
  - given the time it takes to run `get-cross-section.sh`, we try to store the
    most up-to-date version of `xsec.dat` in this repository

### Stream from S3
To stream, we need to make a list of URLs.
- run `generate-s3-list.sh` to generate a list of files
  - running it with no arguments will print the usage and required arguments
  - the file list should appear in `stdout`; pipe the output somewhere, for example:
    - directly to a text file:
      `generate-s3-list.sh S3/.../... > files.txt`
    - add columns for numEvents (0, for all events), cross section (3e-4), and
      Q2min (1), to build a "config" file for an analysis:
      `generate-s3-list.sh S3/.../... 0 3e-4 1 > files.config`
    - use `grep` to remove files that have `"OLD"` in their filename:
      `generate-s3-list.sh S3/.../... 0 3e-4 1 | grep -v OLD > files.config`

### Download from S3
Instead of URLs, we make a list of local files, together with the columns needed to
make a config file
- first, to download the files from S3 to a local directory
  `/my/local/directory`, pipe the output of `generate-s3-list.sh` (see above)
  to `download.sh`, for example:
  - `generate-s3-list.sh S3/.../... | grep -v OLD | download.sh /my/local/directory`
  - this will generate a downloader script `/my/local/directory/get-files.sh`,
    and execute it; this `get-files.sh` script is useful because it tells you
    where the files have been downloaded from, and can be easily executed again
    or editted in case the downloading experience any problems
- next run the script `generate-local-list.sh`, which is very similar to
  `generate-s3-list.sh`, but uses files from the specified local directory; see
  above for details and example usage
  - **Important**: execute this script from the main repository directory
    (`../`), so that the file names are relative to the same directory that the
    macros are designed to be executed from:
    `s3tools/generate-local-list.sh path/to/data`
    - or just specify an absolute path, which is more robust
