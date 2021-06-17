# largex-eic

- set env vars before doing anything
  `source env.sh /path/to/delphes/repository`
  (if you do not specify a path to `delphes` repo, it will use
  a default path)
- use `exeDelphes.sh` to run `delphes` on a specified input file
  - in the script, you may need to change `exeDelphes` to the proper
    executable, e.g., `DelphesHepMC2` or `DelphesHepMC3`
  - run `exeDelphes.sh` with no arguments for usage guide
- build analysis code with `make`
- run analysis code: `analysis.exe [delphes output root file]`
