# largex-eic

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
- build analysis code with `make`
  - it requires a `root` build as well as `delphes`
- run analysis code: `analysis.exe [delphes output root file]`
  - the output will be a `root` file, filled with `TObjArray`s of
    histograms
  - each `TObjArray` can be for a different subset of events, e.g., 
    different minimum y cuts, so that their histograms can be compared
    and divided
