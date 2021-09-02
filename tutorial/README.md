# Tutorials

Here is a collection of tutoial macros. If you are learning this software,
it is recommended to go through these tutorials in order given below.

- **Note**: execute macros from the `largex-eic` top directory, not from
this tutorial dirctory, e.g., `root -b -q tutorial/analysis_template.C`

- you need to generate or obtain ROOT files in the `datarec/` directory

- currently all of these examples focus on fast simulations

- it is recommended to create two macros:
  - analysis macro, to run an `Analysis` class, which will analyze 
    simulated data in an event loop, and generate a variety of output
    data structures
  - postprocessor macro, to process the output from the analysis macro,
    such as drawing plots

## Examples:

1. Template
  - `analysis_template.C`: minimal analysis macro to demonstrate how
    to run `Analysis`
  - there is no postprocessor macro (see other examples); instead, inspect
    the output root file from the analysis macro, to learn what objects
    are stored
