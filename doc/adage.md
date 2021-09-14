# Adage

Analysis in a Directed Acyclic Graph Environment (ADAGE) is a framework designed to integrate user analysis code as part of a generalized multi-dimensional binning data structure. The data structure is a Directed Acyclic Graph (DAG) that has two types of nodes:
- Bin nodes: hold a one-dimensional bin specification, implemented by the `CutDef` class
- Control nodes: hold lambda expressions, which are executable during DAG traversal algorithms

Nodes are implemented by the `Node` class, the DAG is implemented by `DAG`. Nodes are connected by single-direction edges, and any path of nodes through the DAG is handled by the `NodePath` class. Classes that derive from `DAG` allow for the attachment of any data structure to the multidimensional bin scheme; `HistosDAG` is an example which tracks a `Histos` object (a set of histograms) for each multi-dimensional bin.

The figure below shows a DAG for a 4-dimensional binning scheme in (x,y,Q2,z). The small blue circles are control nodes, and the large green rectangles are bin nodes. There is a unique control node at the top, called the root node, and a single control node at the bottom, called a leaf node. Each layer of bin nodes corresponds to one variable, and thus to one dimension. Each layer is fully connected to its neighboring layers; a full connection between two layers is called a "full patch".

![fig1](img/dag1.png)
