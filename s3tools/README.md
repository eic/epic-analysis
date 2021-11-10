# S3 Tools

Scripts to provide access to S3 and interface with this SIDIS analysis
repository

For more details, see [S3 file storage documentation](https://doc.athena-eic.org/en/latest/howto/s3_file_storage.html)

### Quick Start
Here are commands I'm using for testing; documentation below describes what they do
```
# download files for 5x41, 3 different Q2 minima
s3tools/generate-s3-list.sh S3/eictest/ATHENA/RECO/canyonlands-v1.2/DIS/NC/5x41/minQ2=1   | s3tools/download.sh datarec/canyonlands/5x41/minQ2=1
s3tools/generate-s3-list.sh S3/eictest/ATHENA/RECO/canyonlands-v1.2/DIS/NC/5x41/minQ2=10  | s3tools/download.sh datarec/canyonlands/5x41/minQ2=10
s3tools/generate-s3-list.sh S3/eictest/ATHENA/RECO/canyonlands-v1.2/DIS/NC/5x41/minQ2=100 | s3tools/download.sh datarec/canyonlands/5x41/minQ2=100
tree datarec/canyonlands

# build a config file (TODO: need to set $crossSection for each case); Q2 min should be decreasing
s3tools/generate-local-list.sh datarec/canyonlands/5x41/minQ2=100 0 $crossSection 100 | tee    datarec/canyonlands/5x41/files.config
s3tools/generate-local-list.sh datarec/canyonlands/5x41/minQ2=10  0 $crossSection 10  | tee -a datarec/canyonlands/5x41/files.config
s3tools/generate-local-list.sh datarec/canyonlands/5x41/minQ2=1   0 $crossSection 1   | tee -a datarec/canyonlands/5x41/files.config
more datarec/canyonlands/5x41/files.config

# run a macro
root -b -q macro/ci/analysis_dd4hep.C'("datarec/canyonlands/5x41/files.config",5,41)'
```

### Accessing S3 Files
- first, download the [MinIO client](https://docs.min.io/docs/minio-client-complete-guide)
  - if you are using the Singularity or Docker container, it is already installed
  - the main command is `mc`
- next, setup your client to connect to our S3 host (ask someone for credentials):
  - first set env vars `$S3_ACCESS_KEY` and `$S3_SECRET_KEY` with the login and password
    (ask someone if you don't know):
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

# Generating Config Files
Next we need to make a "config file", which consists of the file name, and
additional comments such as Q2min (see [documentation here](../tutorial/README.md)).
Follow the next sections, whether you plan to stream from S3 or download.

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
