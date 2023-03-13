# Tutorials

Here is a collection of tutorial macros. If you are learning this software,
it is recommended to go through these tutorials in the numerical order given below.

- **Note**: Execute macros from the `epic-analysis` top directory, not from
this tutorial directory, e.g., `root -b -q tutorial/analysis_template.C`

## Generate or Obtain Simulation Output Trees 

To run tutorials, you need to generate or obtain ROOT files, from fast or full simulation
- The `datarec/` directory is provided for storage of these files,
  but is not required to use
- The following sub-sections describe how to obtain these files from the common
  storage area
  - This requires access to S3, the common storage area: follow
    [s3tools documentation](../s3tools/README.md) for guidance
- A set of files, together with settings such as beam energy and Q2 ranges, are
  specified by config files; see [doc/example.config](../doc/example.config) for an example
  config file and more details

### Fast Simulation
- To download sample HEPMC files from S3, and run them through Delphes, run:
  ```bash
  s3tools/s3tool.rb -v hepmc.pythia8 -o tutorial.fastsim -c tutorial/delphes.config -e 10x100 -l 4
  ```
  - Run `s3tools/s3tool.rb` (no arguments) for a usage guide of this script
    - Different `hepmc` datasets may be available (control with the `-v` option)
  - Delphes output files will be written to `datarec/tutorial.fastsim`
    (and the HEPMC files will be in `datagen/tutorial.fastsim`)
  - By default, all of the tutorial macros for fast simulations assume the
    `config` file is `tutorial/delphes.config` (as specified by the `-c` option)


### Full Simulation
- The same `s3tools/s3tool.rb` script is used to download full simulation files from S3:
  ```bash
  s3tools/s3tool.rb -e 18x275 -o tutorial.epic -c tutorial/s3files.epic.config -l 4
  ```
  - This is for the latest ePIC data (with the specified beam energy, and
    writes the config file to `tutorial/s3files.epic.config`)
    - By default, the `config` file will be filled with URLs for streaming data
      from S3; if you would rather download the files locally, add the option
      `-m d` to copy them to `datarec/tutorial.epic/`
  - Run `s3tools/s3tool.rb -e 18x275 -v PRODUCTION_VERSION` for a different
    production, such as one from ECCE or ATHENA (run `s3tools/s3tool.rb` with
    no arguments to see available `PRODUCTION_VERSION`s)
  - Similar to the fast simulations, note where the `config` file is produced; the tutorial
    macros require this file as an argument; you can control this with the `-c` option:
    - `s3tools/s3tool.rb -c tutorial/s3files.epic.config` for ePIC macros
    - `s3tools/s3tool.rb -c tutorial/s3files.ecce.config` for ECCE macros
    - `s3tools/s3tool.rb -c tutorial/s3files.athena.config` for ATHENA macros


## Introductory Notes

### Switching between Fast and Full Simulations
- many of these examples focus on fast simulations; to switch between fast and
  full simulations, change the `Analysis`-derived class in the macro:
  - `AnalysisDelphes` for Delphes trees (fast simulations)
  - `AnalysisEpic` for trees from the DD4hep+EICrecon stack (ePIC full simulations)
  - `AnalysisAthena` for trees from the DD4hep+Juggler stack (ATHENA full simulations)
  - `AnalysisEcce` for trees from the Fun4all+EventEvaluator stack (ECCE full simulations)
- note: some extra settings and features differ between these

### Input File Lists (Config Files)
- in the analysis macros, the input files are specified by a list, a "config
  file", with the additional columns
  - see [documentation here](../s3tools/README.md) for the formatting of these
    files, as well as scripts to help generate them
  - for example, this file allows one to combine different Q2 regions together
    using relative weights 


# Tutorials:

Each of these examples has two macros:
  - analysis macro, to run an `Analysis`-derived class, which will analyze 
    simulated data in an event loop, and generate a variety of output
    data structures
  - postprocessor macro, to process the output from the analysis macro,
    such as drawing plots
  - the analysis macro will take some time to run, since it runs
    the event loop; the postprocessor macro is typically fast, since
    it analyzes the resulting data structures


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
    - this example also describes how cuts are defined
  - `postprocess_xqbins_draw.C`: draws a couple sample histograms for
    the given binning scheme

3. Full Simulations (all other tutorials are for fast simulations)
  - `analysis_epic.C`: basically a copy of `analysis_xqbins.C`,
    but shows how to analyze full simulation data; the main difference
    is using `AnalysisEpic` instead of `AnalysisDelphes`
  - `postprocess_epic_draw.C`: clone of `postprocess_xqbins_draw.C`,
    specific for this example
  - see also `analysis_eventEvaluator.C` and `postprocess_eventEvaluator_draw.C`
    for similar full simulation scripts using the `EventEvaluator` output from
    ECCE simulations
  - see also `analysis_athena.C` for the ATHENA version

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
