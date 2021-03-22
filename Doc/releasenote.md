# Release Highlights of PlatEMO 3.1

* Add two multi-objective optimization algorithms CCGDE3 and NSGA-II+ARSBX and one single-objective optimization algorithm OFA. There are currently 153 algorithms in the platform.

* Fix some minor bugs in algorithms and the GUI.

# Release Highlights of PlatEMO 3.0

* 20+ algorithms and 100+ problems for single-objective optimization. There are currently 150 algorithms and 339 problems in the platform, including single-objective optimization, multi-objective optimization, many-objective optimization, combinatorial optimization, large-scale optimization, constrained optimization, multimodal optimization, expensive optimization, sparse optimization, and preference optimization.

* A totally new GUI with more powerful functions, which contains a test module, an application module, and an experiment module.

* A novel filter system based on hybrid labels, which facilitates the selection of suitable algorithms for solving different types of problems.

* More convenient interfaces for solving user-defined problems, where no file needs to be written by users.

* A better visualization of populations, where the true Pareto fronts and feasible regions can be shown in the plots.

# Release Highlights of PlatEMO 2.9

* Add one algorithm for constrained optimization (i.e., CMOEA-MS), one algorithm for large-scale optimization (i.e., DGEA), one algorithm for expensive optimization (i.e., MESMO), and one algorithm for feature selection (i.e., DAEA). There are currently 122 algorithms in the platform.

# Release Highlights of PlatEMO 2.8  

* Add three algorithms for constrained optimization (i.e., CCMO, MOEA/D-DAE, and TiGE-2) and an algorithm for many-objective optimization (i.e., PREA). There are currently 118 algorithms in the platform.
* Fix some minor bugs in the Pareto front sampling methods in LIR-CMOP and MW problems.

# Release Highlights of PlatEMO 2.7  

* Add two algorithms for large-scale optimization, i.e.,  GLMO and LCSA. There are currently 114 algorithms in the platform.
* Add four sparse multi-objective optimization problems, i.e., community detection, instance selection, portfolio optimization, and sparse signal reconstruction.
* Add six constrained DTLZ problems, i.e., DC-DTLZ. There are currently 217 problems in the platform.

# Release Highlights of PlatEMO 2.6  

* Add two algorithms: MOEA/PSL and DWU. There are currently 112 algorithms in the platform.

# Release Highlights of PlatEMO 2.5  
* Add the time-varying ratio error estimation (TREE) test suite, which contains six constrained large-scale problems from real-world applications.
* Fix some minor bugs in algorithms and problems.


# Release Highlights of PlatEMO 2.4  
* Add two algorithms: MSEA and OSP-NSDE. There are currently 110 algorithms in the platform.

# Release Highlights of PlatEMO 2.3  
* Add four algorithms: C-TAEA, ToP, MOEA/D-URAW, and MultiObjectiveEGO. There are currently 108 algorithms in the platform.
* Add the constrained benchmark problems DOC1-9 and MW1-14. There are currently 201 problems in the platform.
* Update the Pareto front sampling methods of DAS-CMOP1-9 and LIR-CMOP1-14: Dynamically sample points on Pareto fronts instead of loading points from files.
* Update the table in experiment module: Ignore NaN values when calculating the mean and standard deviation in each cell of the table.

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
