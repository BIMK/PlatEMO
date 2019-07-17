

# Release Highlights of PlatEMO 2.2  

* Add two algorithms AGE-MOEA and PPS.
* Add the constrained benchmark problems DAS-CMOP1-9 and LIR-CMOP1-14.

# Release Highlights of PlatEMO 2.1
* Add the sparse multi-objective evolutionary algorithm SparseEA.
* Add the sparse multi-objective test suite SMOP1-SMOP8.
* Add four sparse multi-objective optimization problems, i.e., feature  selection, pattern mining, critical node detection, and neural network training.
* Add the diversity metric CPF (i.e., coverage over Pareto front).
* Add the irregular multi-objective test suite IMOP1-IMOP8.

# Release Highlights of PlatEMO 2.0
* __Lighter framework.__ The architecture of PlatEMO is simplified, which leads to lower learning cost and higher efficiency. The result file size is also reduced.  
* __Higher efficiency.__ The runtime of Pareto dominance based algorithms is reduced by using a more efficient non-dominated sorting algorithm. The runtime of decomposition based algorithms is reduced due to the new architecture of PlatEMO. The runtime of hypervolume calculation is reduced by new logic and GPU acceleration. In experimental module, the algorithms can be executed in parallel.  
* __More conveniences.__ The populations obtained during the evolutionary process can be saved in result files. The references of each algorithm, problem, operator, and metric are given in the comments of the function. The codes of GUI are now open source.