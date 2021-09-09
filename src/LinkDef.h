#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class CutDef+;
#pragma link C++ class BinSet+;
#pragma link C++ class Node+;
#pragma link C++ class NodePath+;
#pragma link C++ class DAG+;
#pragma link C++ class HistosDAG+;
#pragma link C++ class HistConfig+;
#pragma link C++ class Histos+;
#pragma link C++ class Kinematics+;
#pragma link C++ class SimpleTree+;
#pragma link C++ class AnalysisDelphes+;
#pragma link C++ class AnalysisDD4hep+;
#pragma link C++ class PostProcessor+;
#pragma link C++ class Weights+;
#pragma link C++ class WeightsUniform+;
#pragma link C++ class WeightsSivers+;
#pragma link C++ class WeightsCollins+;
#pragma link C++ class WeightsProduct+;
#pragma link C++ class WeightsSum+;

// links for DAG streaming: (this clears warnings, but not errors...)
/*
#pragma link C++ class function<void(Node*,NodePath)>+;
#pragma link C++ class _Maybe_unary_or_binary_function<void,Node*,set<Node*>>+;
#pragma link C++ class binary_function<Node*,set<Node*>,void>+;
#pragma link C++ class _Function_base+;
*/


#endif
