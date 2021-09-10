# Tutorials

Here is a collection of tutoial macros. If you are learning this software,
it is recommended to go through these tutorials in order given below.

- **Note**: execute macros from the `largex-eic` top directory, not from
this tutorial dirctory, e.g., `root -b -q tutorial/analysis_template.C`

- you need to generate or obtain ROOT files in the `datarec/` directory

- currently all of these examples focus on fast simulations

- in general it is recommended to use two macros:
  - analysis macro, to run an `Analysis` class, which will analyze 
    simulated data in an event loop, and generate a variety of output
    data structures
  - postprocessor macro, to process the output from the analysis macro,
    such as drawing plots

## Examples:

1. Template
  - `analysis_template.C`: minimal analysis macro to demonstrate how
    to run `Analysis`; no bins are specified
  - there is no postprocessor macro (see other examples); instead, inspect
    the output root file from the analysis macro, to learn what objects
    are stored

2. (x,Q2) Binning
  - `analysis_xqbins.C`: bin the analysis in a few 2D bins of x and Q2
    - there is a `switch` statement to allow the choice of various
      example binning schemes
  - `postprocess_xqbins_draw.C`: draws a couple sample histograms for
    the given binning scheme

3. Average kinematics table
  - `analysis_qbins.C`: bin the analysis in several Q2 bins, for a couple
    pT bins
  - `postprocess_qbins_tables.C`: produce a text file containing tables
    of averages of kinematics variables, where each row is for a Q2 bin;
    one table is produced for each pT bin

4. Full Simulations
  - `analysis_dd4hep_v0.C` and `postprocess_dd4hep_v0.C`
  - all other tutorials are fast simulation; this one shows how to run
    full simulations on some sample files on `S3`
    - note you need to set your `S3` access credentials via environment
      variables `S3_ACCESS_KEY` and `S3_SECRET_KEY`
    - primary difference is the usage of `AnalysisDD4hep` instead of
      `Analysis`

More examples will be added eventually; for now you are encouraged to
look at other existing analysis and postprocessor macros

