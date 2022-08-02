#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class CutDef+;   // TODO: check for dependency on sidis-eic
#pragma link C++ class BinSet+;   // TODO: check for dependency on sidis-eic
#pragma link C++ class Node+;     // TODO: check for dependency on sidis-eic
#pragma link C++ class NodePath+; // TODO: check for dependency on sidis-eic
#pragma link C++ class DAG+;      // DONE: DAG and Adage are independent of sidis-eic

#pragma link C++ class Histos+;
#pragma link C++ class Adage<Histos>+;
#pragma link C++ class HistosDAG+;

#pragma link C++ class HistConfig+;
#pragma link C++ class Hist4D+;
#pragma link C++ class Kinematics+;
#pragma link C++ class SimpleTree+;
#pragma link C++ class Analysis+;
#pragma link C++ class AnalysisDelphes+;
#pragma link C++ class AnalysisAthena+;
#pragma link C++ class AnalysisEcce+;
#pragma link C++ class PostProcessor+;
#pragma link C++ class Weights+;
#pragma link C++ class WeightsUniform+;
#pragma link C++ class WeightsSivers+;
#pragma link C++ class WeightsCollins+;
#pragma link C++ class WeightsProduct+;
#pragma link C++ class WeightsSum+;

#endif
