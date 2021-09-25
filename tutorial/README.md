# Tutorials

Here is a collection of tutoial macros. If you are learning this software,
it is recommended to go through these tutorials in order given below.

- **Note**: execute macros from the `largex-eic` top directory, not from
this tutorial directory, e.g., `root -b -q tutorial/analysis_template.C`

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

3. Full Simulations (all other tutorials are for fast simulations)
  - `analysis_dd4hep.C`: basically a copy of `analysis_xqbins.C`,
    but shows how to analyze full simulation data; the main difference
    is using `AnalysisDD4hep` instead of `AnalysisDelphes`
  - `postprocess_dd4hep_draw.C`: clone of `postprocess_xqbins_draw.C`,
    specific for this example

4. Average kinematics table
  - `analysis_qbins.C`: bin the analysis in several Q2 bins, for a couple
    pT bins
  - `postprocess_qbins_tables.C`: produce a text file containing tables
    of averages of kinematics variables, where each row is for a Q2 bin;
    one table is produced for each pT bin

5. Draw ratio of histograms for different y-minima
  - `analysis_yRatio.C`: bins in 3 y-minima, along with a full-y bin
  - `postprocess_yRatio.C`: produces ratios of all histograms, with
    a y-minimum divided by no y-minimum

6. Test DAG lambda staging
  - `analysis_testDAG.C`: define multi-dimensional binning
  - `postprocess_testDAG.C`: contains various lambdas examples and
    describes how to stage them

7. Conditional Subloops
  - `analysis_coverage.C`: define 4-D binning scheme, including
    extra "full-range" bins
  - `postprocess_coverage.C`: draw certain plots, while restricting
    certain bins to be "full-range", i.e., "integrated over"; in order
    to restrict the execution of certain subloops, conditional control
    functions are used

More examples will be added eventually; for now you are encouraged to
look at other existing analysis and postprocessor macros

