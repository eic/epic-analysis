# Tutorials

Here is a collection of tutorial macros. If you are learning this software,
it is recommended to go through these tutorials in the numerical order given below.

**Note**: Execute macros from the `epic-analysis` top directory, not from
this tutorial directory, for example:
```bash
root -b -q -l tutorial/analysis_template.C
```

## Generate or Obtain Simulation Output Trees 

To run tutorials, you need to generate or obtain ROOT files, from fast or full simulation.
- Follow [the s3tools documentation](../s3tools/README.md) for guidance;
  use the example `s3tools/s3tool.rb` commands within to get some files for these tutorials,
  or write your own `s3tools/s3tool.rb` command to get any other files
- The `datarec/` directory is provided for storage of these files, but is not required to use

## Introductory Notes

### Macros: Analysis and Postprocessing
In general, two macros are used to run `epic-analysis`:
1. Analysis macro, to run an `Analysis`-derived class, which will analyze 
   simulated data in an event loop, and generate a variety of output
   data structures
2. Postprocessor macro, to process the output from the analysis macro,
   such as drawing plots

The analysis macro will take some time to run, since it runs the event loop;
the postprocessor macro is typically fast, since it analyzes the resulting data
structures

### Input File Lists (Config Files)
A set of files, together with settings such as beam energy and Q2 ranges, are
specified by config files; these files are the _input_ to the analysis macros.
See [doc/example.config](../doc/example.config) for an example config file and
more details. The `s3tools/s3tool.rb` script will generate such files for you
automatically.

### Switching between Fast and Full Simulations
- many of these examples focus on fast simulations; to switch between fast and
  full simulations, change the `Analysis`-derived class in the macro:
  - `AnalysisDelphes` for Delphes trees (fast simulations)
  - `AnalysisEpic` for trees from the DD4hep+EICrecon stack (ePIC full simulations)
  - `AnalysisAthena` for trees from the DD4hep+Juggler stack (ATHENA full simulations)
  - `AnalysisEcce` for trees from the Fun4all+EventEvaluator stack (ECCE full simulations)
- note: some extra settings and features differ between these

-----------------------
-----------------------
Tutorial Example Macros
=======================
-----------------------

# 1. Template
- `analysis_template.C`: minimal analysis macro to demonstrate how
  to run `Analysis`; no bins are specified
- there is no postprocessor macro (see other examples); instead, inspect
  the output root file from the analysis macro, to learn what objects
  are stored

Commands:
```bash
root -b -q -l tutorial/analysis_template.C   # run the analysis macro
root out/tutorial.template.root              # open the ROOT file
```

Once in the ROOT file, try:
```cpp
new TBrowser               // open a browser
tree->Draw("Q2")           // draw a Q2 distribution (not weighted!)
tree->Draw("Q2","Weight")  // draw a Q2 distribution, with the correct weights
```
Explore the ROOT file, view the histograms in the `histArr*` TDirectories.

Try to run on ePIC full simulations; change the macro to run with `AnalysisEpic`,
then call
```bash
root -b -q -l tutorial/analysis_template.C'("tutorial/epic.config")'
root out/tutorial.template.root
```

# 2. (x,Q2) Binning
- `analysis_xqbins.C`: bin the analysis in a few 2D bins of x and Q2
  - there is a `switch` statement to allow the choice of various
    example binning schemes
  - this example also describes how cuts are defined
- `postprocess_xqbins_draw.C`: draws a couple sample histograms for
  the given binning scheme

Commands:
```bash
root -b -q -l tutorial/analysis_xqbins.C           # run the analysis macro
root -b -q -l tutorial/postprocess_xqbins_draw.C   # run the post-processing macro
root out/tutorial.xqbins.canvas.root               # open the resulting ROOT file
```
View the images in `out/tutorial.xqbins.images/`, for example:
```bash
display out/tutorial.xqbins.images/*.png   # press 'space' or 'backspace' to change images
```

# 3. Full Simulations
All other tutorials are for fast simulations (by default); this one is an example for full simulations
- `analysis_epic.C`: basically a copy of `analysis_xqbins.C`,
  but shows how to analyze full simulation data; the main difference
  is using `AnalysisEpic` instead of `AnalysisDelphes`
- `postprocess_epic_draw.C`: clone of `postprocess_xqbins_draw.C`,
  specific for this example
- see also `analysis_eventEvaluator.C` and `postprocess_eventEvaluator_draw.C`
  for similar full simulation scripts using the `EventEvaluator` output from
  ECCE simulations
- see also `analysis_athena.C` for the ATHENA version

# 4. Average kinematics table
- `analysis_qbins.C`: bin the analysis in several Q2 bins, for a couple
  pT bins
- `postprocess_qbins_tables.C`: produce a text file containing tables
  of averages of kinematics variables, where each row is for a Q2 bin;
  one table is produced for each pT bin

# 5. Draw ratio of histograms for different y-minima
- `analysis_yRatio.C`: bins in 3 y-minima, along with a full-y bin
- `postprocess_yRatio.C`: produces ratios of all histograms, with
  a y-minimum divided by no y-minimum

# 6. Test DAG lambda staging
- `analysis_testDAG.C`: define multi-dimensional binning
- `postprocess_testDAG.C`: contains various lambdas examples and
  describes how to stage them

# 7. Conditional Subloops
- `analysis_coverage.C`: define 4-D binning scheme, including
  extra "full-range" bins
- `postprocess_coverage.C`: draw certain plots, while restricting
  certain bins to be "full-range", i.e., "integrated over"; in order
  to restrict the execution of certain subloops, conditional control
  functions are used

# 8. PODIO Examples
Additional example macros demonstrating PODIO usage:
- `podio_frame_reader.C`
  - simple example showing how to read a ROOT file with PODIO
  - `Analysis`-derived classes which use PODIO are essentially based
    of this procedure
- `podio_print_collection_datatypes.C`
  - dump the PODIO Collection type for all Collections contained in
    a ROOT file
  - it is very useful to know what datatypes and data models are available
