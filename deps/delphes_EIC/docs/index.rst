.. Documentation master file

Delphes EIC
===========

.. toctree::
   :maxdepth: 2
   :caption: Getting started:

   install.md

   
.. toctree::
   :maxdepth: 2
   :caption: Links:

   GitHub <https://github.com/eic/delphes_EIC>
   Delphes <https://github.com/delphes/delphes>
   
   
.. Add main page text here


Welcome to the **Delphes EIC** project. This code aims to deliver a fast-simulation model of baseline or specific detectors for studies to support the Electron-Ion Collider (EIC) project. This code has been used in a few publications and papers and parts of it are available for citation through Zenodo:

- Charm jets as a probe for strangeness at the future Electron-Ion Collider. Miguel Arratia (UC Riverside amd Jefferson Lab), Yulia Furletova (Jefferson Lab), T. J. Hobbs (Southern Methodist U. and Jefferson Lab), Frederick Olness (Southern Methodist U.), Stephen J. Sekula (Southern Methodist U.). https://arxiv.org/abs/2006.12520 [hep-ph]
- A Delphes card for the EIC yellow-report detector. Miguel Arratia (UC Riverside amd Jefferson Lab) and Stephen J. Sekula (Southern Methodist U.). https://arxiv.org/abs/2103.06886 [physics.ins-det].  DOI: https://doi.org/10.5281/zenodo.4592887.
- Science Requirements and Detector Concepts for the Electron-Ion Collider: EIC Yellow Report.R. Abdul Khalek, A. Accardi, J. Adam, D. Adamiak, W. Akers et al. e-Print: https://arxiv.org/abs/2103.05419 [physics.ins-det]


.. image:: https://github.com/eic/delphes_EIC/raw/master/images/EICDetector_3D_CCDIS_CharmJet.png

.. image:: https://github.com/eic/delphes_EIC/raw/master/images/EICDetector_3D_CCDIS_CharmJet_DisplacedVtx.png

====
What's New?
====

* EIC PID code has been used to create IdentificationMaps for the mRICH, barrelDIRC, and dualRICH. No more external code needed!

====
Installation Instructions
====

#. OPTIONAL (but recommended): Install LHAPDF
#. Install PYTHIA8.3 (if you installed LHAPDF, build with LHAPDF support)
#. Install Delphes3 following: https://github.com/delphes/delphes


====
EIC Yellow-Report Detector Models
====

The current model we recommend is ```delphes_card_allsilicon_3T.tcl```. Is is based on detector studies from the EIC Yellow Report (https://arxiv.org/abs/2103.05419). This model incorporates all-silicon tracking, as well as an EMCAL and HCAL. PID system responses are provided by efficiency maps based on EIC PID code (https://gitlab.com/preghenella/pid/). 

- Magnetic field: 1.5 T or 3.0 T
- Solenoid length: 2.0 m
- Tracker radius: 80 cm
- An HCAL and and EMCAL
- Particle ID systems: based on the mRICH, barrel DIRC, and dual RICH concepts articulated by the EIC community; implemented using IdentificationMaps

We currently simulate DIS using Pythia8 within Delphes. The command file (ending in `.cmnd`) shown here and available in the project is suitable for DIS at EIC. 

====
Using the project
====

Generating Events
----

We recommend PYTHIA 8.3  for event generation, and the latest version of Delphes builds and links against 8.3 without incident, allowing for seamless generation and then detector response modeling all in one step using the ```DelphesPythia8``` binary:

.. code:: bash

   DelphesPythia8 delphes_card_allsilicon_3T.tcl pythia8cards/CC_DIS.cmnd out.root

For examples of other detector configurations, analysis code, etc. see the main Delphes project site. 

Visualizing Events 
----

Run the Delphes visualization script on your ROOT file, using the input Delphes card (TCL) file to help it structure the detector in the display:

.. code:: bash

   root -l examples/EventDisplay.C'("delphes_card_allsilicon_3T.tcl","out.root")'
 
The two examples shown here are for neutral-current and charged-current event for beam energies of 10 GeV electron on 100 GeV proton (63 GeV center-of-mass energy). 

EIC Collider Variations
----

Beam energy recommended benchmarking points are (the order is hadron on lepton):

* 275 on 18 GeV
* 100 on 10 GeV
* 100 on 5 GeV
* 41 on 5 GeV


The "SimpleAnalysis" framework
----

See the dedicated SimpleAnalysis_ Documentation.

.. _SimpleAnalysis: SimpleAnalysis/README.md
