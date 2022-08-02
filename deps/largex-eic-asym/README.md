# largex-eic-asym

## Building

- clone with submodules 
  `git clone --recurse-submodules <this-repo's-clone-address>`
- build BruFit
  - `source env.sh`
  - `pushd deps/brufit`
  - `mkdir build`
  - `cd build`
  - `cmake ../`
  - `make install` (probably will have warnings, since on dev)
  - `cmake ../` (again, for pcm file)
  - `make install`
  - `popd`
- build Largex-eic-asym
  - `make` (make sure `env.sh` has been sourced)

## Usage

- `source env.sh`
- put a `root` file with a `SimpleTree` TTree in `data/` (or anywhere you want)
- run `root -b -q indexTree.C'(...)'` to index the tree
  - this will make a copy with the tree, with the file extension `.root` modified
    to `.idx.root`
  - the new tree will have a branch `Idx`, a unique event index needed for the
    weights implementation
- run `brufit -b -q pullWeights.C'(...)'` to read the specified weights branch from the
  tree, and produce an output file with the `Weights` data structure
  - **important**: run this with `brufit`, not `root`
- edit `brufit.C` and run with `brufit -b -q brufit.C'(...)'`
  - specify a directory name for output files
  - see in particular the `LoadDataSets` call, to specify file and weights branch names;
    if you do not want to use weights, do not specify a weights file
  - use `quantiles.C` to help define binning, if you want quantile binning; 
  - specify which modulations you want to simultaneously fit for; see `Modulation`
    class for more information
  - there are three minimizers available:
    - `minuit` uses Minuit MIGRAD to search for the maximum likelihood
    - `mcmc` uses Metropolis-Hastings Markov Chain Monte Carlo to sample the posterior
      near the likelihood
    - `mcmcthencov` uses the covariance from one MCMC as guidance for a second
  - uses PROOF for multithreading (one bin = one job); see `BruAsymmetry::Fit`
    to disable PROOF and run single-threaded
  - several warnings may be printed, related to "explicitly defaulted constructor"; this is
    a Brufit issue and needs to eventually be fixed
- run `errorCheck.sh` to grep for errors in PROOF log files
  - if you don't use PROOF, check `stdout` for errors
  - don't get confused by uncertainty reports; hopefully the implemented `grep` filters them out
  - some errors are benign
  - make sure `$PROOF_LOG` is correctly set for your system in `env.sh`
- run `draw.C` (this was likely called automatically by `brufit.C`)
  - a new root file will appear in the output directory, containing the fit
    results, along with several other plots
