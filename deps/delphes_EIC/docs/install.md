## Installation


### EJPM install

ejpm is EIC software centric package/build manager. It is the designed
to be the default tool to build G4E on users machine, as it helps with:

- building dependent packages (with "right" compilation flags)
- setup environment variables to run everything
- help with what packages to install by system packet manager
- rebuild/remove/clean existing packages

> Still, ejpm is not a requirement.

First, install ejpm itself:

```bash
   pip install --user ejpm
```
  


(If you have certificate problems (JLab issue), don't have pip, or have other problems,
[here is the detailed documentation](https://gitlab.com/eic/ejpm)

Install delphes and possible other dependencies packets (lhapdf, pythia8, fastjet)


```bash

   # 1. System prerequesties (the most prerequesties are for CERN ROOT)
   ejpm req centos          # get list of required OS packets. Use `ubuntu` on debian  
   sudo yum install ...     # install watever 'ejpm req' tells you

   # 2. Where to install
   ejpm --top-dir=<where-to>   # Directory where packets will be installed

   # 4. Install lhapdf, pythia9 and delphes
   ejpm install delphes

   # 5.  Source environment
   source ~/.local/share/ejpm/env.sh  # Use *.csh file for tcsh

```




You have ROOT and Geant4 and don't want EJPM to build them?  
(Use your installations of ROOT and Geant4)

```bash

   # Before running 'ejpm install g4e'
   ejpm set root `$ROOTSYS`    # Path to ROOT installation
   ejpm set geant <path>       # Path to Geant4 installation   
```

> (!) If you use your version of ROOT, all packages depending on ROOT should be
> installed with the same C++ standard flags as root. So it it was C++11 or C++17, it should be used
> everywhere. To set it in ejpm  
> ```ejpm config global cxx_standard=17```
>

Hint (!). Run ejpm to overview all installed packets, environment and status by 'ejpm' command

Here is the sample output:

```

   > ejpm

   EJPM v0.01.19
   top dir :
      /eic
   state db :
      ~/.local/share/ejpm/db.json  (users are encouraged to inspect/edit it)
   env files :
      ~/.local/share/ejpm/env.sh
      ~/.local/share/ejpm/env.csh

   INSTALLED PACKETS: (*-active):
   vgm:
      * /eic/vgm/vgm-v4-5 (owned) 
   root:
      * /eic/root/root-v6-16-00 (owned)
   geant:
         /eic/geant/geant-v10.5.0 (owned)
      * /eic/geant4-10.6-betta
   hepmc:
      * /eic/hepmc/hepmc-HEPMC_02_06_09 (owned) 
   g4e:
      * /eic/g4e/g4e-dev (owned)
```


### Manual installation

Below, the environment variable ${INSTALLDIR} refers to some folder where you are putting all your EIC fast simulation code (e.g. export INSTALLDIR=${HOME}/EIC/).

1. Install PYTHIA8,
   * http://home.thep.lu.se/~torbjorn/Pythia.html,
   * Download the tarball and unpack it. ,
   * Configure it for local installation in your work area, make, install
   * If you intend to use PDFs besides the default PYTHIA PDF, you need to install LHAPDF6. The configuration instructions below assume you have done this already and want to build PYTHIA with support for LHAPDF
   ```
   ./configure --prefix=${INSTALLDIR}/EIC/ --with-lhapdf6=${INSTALLDIR}/EIC/
   make -j
   make install
   ```
   * Make sure the work area binary directory is in your PATH: ```PATH=${INSTALLDIR}/EIC/bin:${PATH}```,

1. Install Delphes,
   * https://github.com/delphes/delphes,
   * Clone the project and make sure you are on the master branch,
   * Make sure ROOT is available in your path, e.g. ```lsetup \"root 6.18.04-x86_64-centos7-gcc8-opt\"```,
   * Compile with PYTHIA8: ```HAS_PYTHIA8=true PYTHIA8=${INSTALLDIR}/EIC ./configure --prefix=${INSTALLDIR}/EIC/```,
   * Build, install 
   ```
   make -j
   make install
   ```

1. Get the Delphes/EIC code for simulation and analysis of a detector baseline/configuration.,
   * https://github.com/eic/delphes_EIC,
   * Clone the repository locally,
   * Follow the instructions to run the example and generate a ROOT file.


