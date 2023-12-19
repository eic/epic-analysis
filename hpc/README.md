High Performance Computing (HPC) Support
========================================

To run jobs multi-threaded or on a computing cluster, `epic-analysis` can be configured to run one
job per simulation `ROOT` file. The scripts in this directory support this procedure.

Run any script without any arguments for more documentation.

**NOTE**: Running `epic-analysis` using these `hpc` tools is not as well tested as running
single-threaded; please check everything carefully, especially Q2 weights. Report any
issues; you are welcome to contribute your own scripts to support your preferred computing cluster.
It is highly recommended to test jobs with small samples, before launching a full-scale analysis
on all available data.

# Pipeline Automation

The `hpc` toolkit has a built-in macro for streamlining analysis across many campaigns, Q2 ranges, energy configurations, and even detector setups. The pipeline's aim is to automate the following steps entirely on Jefferson Lab's slurm system:

1. Creation of the main s3 `.config` file (typically stored in `datarec/`)
2. Calculation of the number of events stored within each `s3` file's TTree (used for calculating event-by-event weights). These are also cached in `hpc/nevents_databases` for faster computation of future pipelines.
3. Splitting of the main s3 `.config` file into batches (for parallel computing).
4. Execution of the analysis macro for each batched `.config` file
5. Merging of the output analysis `.root` files into a single `analysis.root` file

The script that handles the pipeline is `run-local-slurm-pipeline.rb`. The user should edit this script with the desired configurations. These include the campaigns, the energies of interest within those campaigns, the detector configuration, the number of files from `s3` to analyze (per Q2 range) and the number of root files which are analyzed per slurm job. By default, several of these parameters will trip the error handler until the user sets them accordingly.

We note that the calculation of the `nevents` for each s3 `TTree`, albeit time-consuming, is very important for our parallel computing needs. This is because the event-by-event Q2weights depend on how many total events are simulated for each Q2 range. Since we are batching the main s3 `.config` into smaller chunks, this information is lost unless we calculate the number of events before running the analysis. These event counts are then used to set manual weights in the batched `.config` files.

To run the pipeline:

```
hpc/run-local-slurm-pipeline.rb
```

Optionally, you can use the `--overwrite` flag to skip the query to delete pre-existing project files.

There are several known issues to be aware of pertaining to memory usage. If `NROOT_FILES_PER_JOB` is too large, then the per job memory allocation listed in `run-local-slurm.rb` may be too small to create the ROOT TTree's from the analysis macro. Additionally, the merging of all ROOT TFile's into one may run out of memory. This would be limited by the memory allocation listed in the pipeline job created by `run-local-slurm-pipeline.rb`. It is set to `4000mb` now which is reasonable and should not run out of memory.

## 1. Preparation
```bash
hpc/prepare.rb
```
A typical `config` file will list several `ROOT` files; run `hpc/prepare.rb` to split a `config`
file into one `config` file per `ROOT` file; these `config` files can each be fed to an analysis
macro. Total yields per Q2 bin are automatically obtained and stored in all `config` files, to make
sure the resulting Q2 weights are correct for the combined set of files.

Alternatively, one can split the starting `config` file into multiple `config` files where the user specifies (as a third argument) the number of `ROOT` files per `config`. To do so, one would utilize the following script.

```bash
hpc/prepare-multi-roots.rb
```




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

### RCF Cluster
- TODO
- Can we avoid using S3?

### JLab Cluster
- TODO - need Slurm config generator
- Can we avoid using S3?

## 3. Merge Output Files
```bash
hpc/merge.rb
```
After successfully running jobs, combine the resulting output `ROOT` files; this is basically `hadd`
but with some handlers for our custom classes `Histos`, `BinSet`, etc. The resulting combined file
can then be used in downstream post-processing macros or user analysis.
