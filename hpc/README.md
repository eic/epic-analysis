High Performance Computing (HPC) Support
========================================

To run jobs multi-threaded or on a computing cluster, `epic-analysis` can be configured to run one
job per simulation `ROOT` file. The scripts in this directory support this procedure.

Run any script without any arguments for more documentation.

Procedure
=========

## 1. Preparation
```bash
hpc/prepare.rb
```
A typical `config` file will list several `ROOT` files; run `hpc/prepare.rb` to split a `config`
file into one `config` file per `ROOT` file; these `config` files can each be fed to an analysis
macro. Total yields per Q2 bin are automatically obtained and stored in all `config` files, to make
sure the resulting Q2 weights are correct for the combined set of files.

## 2. Run Jobs
This step depends on where you want to run jobs. In general, output ROOT files will be written
to a user-specified subdirectory of `out/`.

### Local Condor
```bash
hpc/run-local-condor.rb
```
If you have a local `condor` service, use this script to prepare a `condor` configuration script.
The user must then run `condor_submit` from outside of `eic-shell` (individual jobs will be run in `eic-shell`).
Log files will be written to `hpc/log/`.

## 3. Merge Output Files
```bash
hpc/merge.rb
```
After successfully running jobs, combine the resulting output `ROOT` files; this is basically `hadd`
but with some handlers for our custom classes `Histos`, `BinSet`, etc. The resulting combined file
can then be used in downstream post-processing macros or user analysis.
