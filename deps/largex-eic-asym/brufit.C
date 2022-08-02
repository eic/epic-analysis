R__LOAD_LIBRARY(LargexAsym)
#include "BruAsymmetry.h"

/* arguments
 * - bruDir: name of output directory; several files will be written
 * - minimizer: optimizer algorithm, choose one of the following:
 *   - minuit: maximum likelihood using Minuit MIGRAD
 *   - mcmc: sample posterior using Metropolis-Hastings walk with
 *     gaussian proposal; set hyperparameters at the end of this macro
 *   - mcmcthencov: first run mcmc as above, then use the resulting
 *     covariance matrix in the gaussian proposal function of a second,
 *     subsequent mcmc walk
 */
// - IMPORTANT: execute with `brufit brufit.C'(....)'` (you can do `brufit -b -q`, if you want)

void brufit(
    TString dataTree="../out/asym.idx.root",
    TString bruDir="bruspin",
    TString minimizer="minuit"
    )
{

  // instantiate brufit
  BruAsymmetry * B = new BruAsymmetry(bruDir,minimizer);

  // build modulations
  // - see `Modulation::Initialize` for the full set of available modulations
  B->AddNumerMod(new Modulation(0,Modulation::kUT)); // sin(phiH-phiS)
  B->AddNumerMod(new Modulation(1,Modulation::kUT)); // sin(phiH+phiS)
  B->AddNumerMod(new Modulation(2,Modulation::kUT));
  B->AddNumerMod(new Modulation(3,Modulation::kUT));
  B->AddNumerMod(new Modulation(4,Modulation::kUT));

  // build full PDF
  B->BuildPDF();

  // set binning scheme (see quantiles.C)
  // - bin boundaries (size = number_of_bins + 1)
  Double_t Xbins[7] = { 0.050, 0.061, 0.076, 0.094, 0.119, 0.172, 0.820 };
  // - bin schemes (2nd argument = number_of_bins)
  B->Bin("X",6,Xbins);
  // - check binning definition
  B->PrintBinScheme();

  // load SimpleTrees
  B->LoadDataSets(
      dataTree, // data to fit (must be an indexed tree, if using weights (see indexTree.C))
      "", // unpolarized MC data for likelihood normalization approximation (not used if unspecified)
      "data/Tweights.root", // Tweights file (from pullWeights.C)
      "Weight", // weight branch name
      "tree" // tree name (SimpleTree)
      );

  // optimizer hyperparameters (Markov Chain); not used if minimizer==minuit
  // - MCMC settings
  B->MCMC_iter        = 5000; // number of samples
  B->MCMC_burnin      = 0.1 * ((Double_t)B->MCMC_iter); // number to burn
  B->MCMC_norm        = 1.0 / 0.03; // ~ 1/stepsize
  // - 2nd MCMC settings (if using MCMCthenCov algorithm)
  B->MCMC_cov_iter    = 5000; // number of samples
  B->MCMC_cov_burnin  = 0.1 * ((Double_t)B->MCMC_iter); // number to burn
  B->MCMC_cov_norm    = 1.0 / 0.03; // ~ 1/stepSize

  // perform fit; to decide whether to run single or multi-threaded (with PROOF),
  // see `src/BruAsymmetry.cpp` calls to `HS::FIT::PROCESS::{Here,Proof}::Go`
  B->Fit();

  // print MCMC acceptance rates
  gSystem->RedirectOutput(B->GetLogName());
  gROOT->ProcessLine(".! mcmcAcceptanceRate.sh");
  gSystem->RedirectOutput(0);

  // draw
  TString cmd = Form(".x draw.C(\"%s\",\"%s\")",bruDir.Data(),minimizer.Data());
  printf("\nEXECUTE: %s\n\n",cmd.Data());
  gROOT->ProcessLine(cmd);
};
