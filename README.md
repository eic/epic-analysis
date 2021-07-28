# Largex-eic

# Delphes
- this repository provides a simple wrapper for `delphes`
- first make sure you have a build of `delphes` somewhere, preferably
  in a separate directory
- set env vars before doing anything, so this repo knows where your
  `delphes` build is: `source env.sh /path/to/delphes/repository`
  - if you do not specify a path to `delphes` repo, it will use
  a default path; you can edit this default path for your own 
  convenience
  - it will also symlink `delphes` external code, so analysis macros
    will not complain
- use `exeDelphes.sh` to run `delphes` on a specified input file
  - in the script, you may need to change `exeDelphes` to the proper
    executable, e.g., `DelphesHepMC2` or `DelphesHepMC3`, depending
    on the format of your generator input
  - run `exeDelphes.sh` with no arguments for usage guide
  - this script is a convenience, for automatically handling filenames
    and streaming through `gunzip`, if reading a `.gz` file


# Analysis
- build analysis code with `make` (env vars must be set first (see above))
  - it requires a `root` build as well as `delphes`
  - all classes are stored in the `src/` directory
- the `Analysis` class is the main class that performs the analysis; you 
  can control it from a macro
  - example macros will eventually be added; for now you can assume any macro
    named `analysis_*.C` is an example
  - the analysis macro must the following
    - instantiate `Analysis` (with file names, beam energies, etc.)
    - set up bin schemes and bins (arbitrary specification)
    - set up anything else
    - run analysis
  - the output will be a `root` file, filled with `TObjArray`s of
    histograms
    - each `TObjArray` can be for a different subset of events (bin), e.g.,
      different minimum y cuts, so that their histograms can be compared and
      divided
    - the `Histos` class contains the histograms, and instances of `Histos`
      will also be streamed to `root` files

# PostProcess
- results processing is handled by the `PostProcessor` class, which does tasks
  such as printing tables of average values, and drawing ratios of histograms
  - this class is steered by postprocessor macros
  - example macros will eventually be added; for now you can assume any macro
    named `postprocessor_*.C` is an example
