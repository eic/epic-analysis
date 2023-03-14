# S3 Tools

Scripts to provide access to data on S3 and interface with this SIDIS analysis repository

NOTE: before anything, be sure to set the environment:
```bash
source environ.sh
```

## Setup
First you must obtain access to S3. Set the environment variables `$S3_ACCESS_KEY` and `$S3_SECRET_KEY` with the login and password beforehand:
```bash
export S3_ACCESS_KEY=*****
export S3_SECRET_KEY=*****
```
  - It is useful to add these variables to your shell configuration (such as `~/.bashrc`),
    but do not expose the username and password to the outside world
  - For more details on S3, see [S3 file storage documentation](https://doc.athena-eic.org/en/latest/howto/s3_file_storage.html)

Then add our S3 endpoint:
```bash
s3tools/add-host.sh
```
  - This only needs to be done once, since it will write your configuration
    locally (most likely to `~/.mc/config.json`)

## s3tools Usage
Now you can run our S3 automation script, `s3tools/s3tool.rb`; run without any arguments to print the usage guide:
```bash
s3tools/s3tool.rb
```
- This script can:
  - Download files from S3, given your choice of campaign production version,
    beam energy, detector configuration, radiative corrections, and more
  - Automatically generate a `config` file, with file names, Q2 minima, and
    cross sections; this is used as input to the analysis macros (see
    [doc/example.config](../doc/example.config) for a sample config file)
  - Alternative to downloading files from S3, you can generate a `config` file
    for streaming from S3
- This script supports both fast and full simulations from ePIC, ECCE, and ATHENA
- For fast simulations, event-generated `hepmc` files are obtained and passed through
  `Delphes` locally
- This script uses [MinIO client](https://min.io/docs/minio/linux/reference/minio-mc.html)
  (included in `eic-shell`).

### Example `s3tool.rb` Commands
For full simulation files from S3, run:
```bash
s3tools/s3tool.rb -o tutorial.epic -c tutorial/epic.config -e 18x275 -l 3
```
  - By default, the `config` file will be filled with S3 URLs for streaming data
    from S3; if you would rather download the files locally, add the option
    `-m d` to download them to `datarec/tutorial.epic/`
  - This is for the latest ePIC data, with the specified beam energy, and
    writes the config file to `tutorial/epic.config`
    - Run `s3tools/s3tool.rb -v PRODUCTION_VERSION` for a different
      production, such as one from ECCE or ATHENA (run `s3tools/s3tool.rb` with no
      arguments to see available `PRODUCTION_VERSION`s).

For fast simulation, download sample HEPMC files from S3 and run them through Delphes using:
```bash
s3tools/s3tool.rb -v hepmc.pythia8 -o tutorial.fastsim -c tutorial/delphes.config -e 10x100 -l 3
```
  - **Note:** Delphes must be fully compiled before doing this
  - Delphes output files will be written to `datarec/tutorial.fastsim`,
    and the HEPMC files will be stored in `datagen/tutorial.fastsim`
  - The `config` file will be written to `tutorial/delphes.config`


## Additional Notes

### MinIO Client Usage
If you would rather do things yourself, you can use MinIO client directly using the
`mc` command to browse data on S3
- Example S3 navigation commands:
  - top-level ePIC directory list: `mc ls S3/eictest/EPIC`
  - show directory tree: `mc tree S3/eictest/EPIC/RECO/`
  - download a file: `mc cp /S3/path/to/some/file path/to/local/target/directory/`
  - more documentation: `mc -h`
- Once you have some files, you will need to write your own `config` file
  - Follow [doc/example.config](../doc/example.config) as a template
  - Cross sections are available in [`datarec/xsec/xsec.dat`](../datarec/xsec/xsec.dat)

### Delphes Automation
For convenience, `s3tools/src/loop_run_delphes.sh` can run Delphes on a list of files
in a list of directories
- Takes a list of directories as the arguments; they must be in the `datagen/`
  directory, since the Delphes output files will be to the same path, but with
  `datagen/` replaced with `datarec/`
- Runs multi-threaded: one thread per directory

### Cross Sections
In case you want to update the cross section table `xsec.dat`:
- the script `s3tools/src/get-cross-section.sh` will read the cross section from
  `GenCrossSection` in a `hepmc` file; use `s3tools/src/get-cross-section-ALL.sh` to
  automate running `get-cross-section.sh` over all `hepmc` files in a specific
  directory tree; this will populate `datarec/xsec/*.xsec` files, one for each
  `hepmc` file
- next use `s3tools/src/tabulate-cross-section.py` to read the tree of `.xsec` files into
  a table of cross sections, output to `datarec/xsec/xsec.dat`
- given the time it takes to run `get-cross-section.sh`, we try to store the
  most up-to-date version of `xsec.dat` in this repository
