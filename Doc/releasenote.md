# Release Highlights of PlatEMO 4.1

* Automated function creation is supported. Users can input a dataset as an objective function or constraint function when solving user-defined problems, where a function will be automatically fitted according to the dataset.

* Add two large-scale multi-objective evolutionary algorithms FLEA and LERD, one expensive multi-objective optimization algorithm SMOA, and one constrained multi-objective evolutionary algorithm C3M. There are currently 220 algorithms in the platform.

* Add 16 constrained multi-objective benchmark problems ZXH_CF1-ZXH_CF16. There are currently 448 problems in the platform.

# Release Highlights of PlatEMO 4.0

* Dynamic optimization, multitasking optimization, bilevel optimization, and robust optimization are now supported in PlatEMO.

* Hybrid encoding is now supported in PlatEMO, where a problem can include real variables, integral variables, label variables, binary variables, and permutation variables simultaneously.

* Maximum runtime is provided as a new termination criterion, which can be set instead of maximum number of function evaluations.

* More algorithms and problems for single-objective optimization, multi-objective optimization, constrained optimization, sparse optimization, expensive optimization, multimodal optimization, dynamic optimization, multitasking optimization, bilevel optimization, and robust optimization. There are currently 216 algorithms and 432 problems in the platform.

* More efficient and powerful GUI, where the execution of algorithms in the test module and application module is highly accelerated.

* More performance metrics for different types of optimization problems, and the metrics are also tagged with labels. Different metrics will be shown in the dropdown lists when selecting different labels in the GUI.

* Gradient based search is now supported in PlatEMO, where users can define gradient functions to accelerate the convergence via mathematical programming algorithms and gradient assisted evolutionary algorithms.

# Release Highlights of PlatEMO 3.5

* Enhance the application module, where users can define problems and save results more easily.

* Add three decomposition based multi-objective evolutionary algorithms MOEA/D-DCWV, MOEA/D-PFE, and MOEA/D-VOV and a surrogate-assisted multi-objective evolutionary algorithm MCEA/D. There are currently 180 algorithms in the platform.

# Release Highlights of PlatEMO 3.4

* Remake the application module, a more powerful and friendly interface enables users to define problems more easily. The defined problems can also be saved into files and solved in other modules.

* Add two multi-objective evolutionary algorithms MOEA/D-DYTS and MOEA/D-UR, two surrogate-assisted multi-objective evolutionary algorithms PB-NSGA-III and PB-RVEA, a constrained multi-objective evolutionary algorithm DSPCMDE, three large-scale multi-objective evolutionary algorithms FDV, IM-MOEA/D, and LMOEA-DS, two sparse multi-objective evolutionary algorithm SLMEA and SparseEA2, and four single-objective mathematical programming methods Adam, Nelder-Mead, RMSProp, and SD. There are currently 176 algorithms in the platform.

# Release Highlights of PlatEMO 3.3

* Add four multi-objective evolutionary algorithms DEA-GNG, ICMA, PeEA, and RVEA-iGNG. There are currently 162 algorithms in the platform.

* Add five constrained multi-objective optimization problems FCP1-FCP5 and a sparse multi-objective optimization problem Sparse_KP. There are currently 345 problems in the platform.

# Release Highlights of PlatEMO 3.2

* Add four surrogate-assisted multi-objective evolutionary algorithms AB-SAEA, EDN-ARMOEA, HeE-MOEA, KTA2, and a constrained multi-objective evolutionary algorithm c-DPEA. There are currently 158 algorithms in the platform.

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
