 
  Copyright (C) 2018 Heiner Zille
 
  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

  Author of this Code: 
   Heiner Zille <heiner.zille@ovgu.de>

  This file belongs to the following publications:

  1) Heiner Zille and Sanaz Mostaghim
     "Comparison Study of Large-scale Optimisation Techniques on the LSMOP Benchmark Functions"  
     IEEE Symposium Series on Computational Intelligence (SSCI), IEEE, Honolulu, Hawaii, November 2017
     https://ieeexplore.ieee.org/document/8280974 
 
  2) Heiner Zille, Hisao Ishibuchi, Sanaz Mostaghim and Yusuke Nojima
     "A Framework for Large-scale Multi-objective Optimization based on Problem Transformation"
     IEEE Transactions on Evolutionary Computation, Vol. 22, Issue 2, pp. 260-275, April 2018.
     http://ieeexplore.ieee.org/document/7929324
  
  3) Heiner Zille, Hisao Ishibuchi, Sanaz Mostaghim and Yusuke Nojima
     "Weighted Optimization Framework for Large-scale Mullti-objective Optimization"
     Genetic and Evolutionary Computation Conference (GECCO), ACM, Denver, USA, July 2016
     http://dl.acm.org/citation.cfm?id=2908979


  The source code was programmed to work with the PlatEMO Framework, version 1.5 for Matlab [1]. 
  
  To use this program, simply copy the complete "WOFSMPSO"-folder and its contained files 
  into the "Algorithms"-folder of the PlatEMO framework. WOFSMPSO should then be usable as any other 
  algorithm inside the framework. 
  
  File list:     
    - WOFSMPSO/gpl-3.0.txt
    - WOFSMPSO/README.txt
    
    - WOFSMPSO/WOF_createGroups.m
    - WOFSMPSO/WOF_CrowdingDistance.m
    - WOFSMPSO/WOF_EnvironmentalSelection.m
    - WOFSMPSO/WOF_optimiseBySMPSO.m
    - WOFSMPSO/WOF_selectxPrimes.m
    - WOFSMPSO/WOF_SMPSO_operator.m
    - WOFSMPSO/WOF_transformationFunction.m
    - WOFSMPSO/WOF_WeightIndividual.m
    - WOFSMPSO/WOFSMPSO.m
    
  Date of uploading these files: 12.10.2018 
  Last Update: 12.10.2018

  You can reach me via email at <heiner.zille@ovgu.de> or visit my webpage at: 
  http://www.ci.ovgu.de/Team/Heiner+Zille.html
  
  
  [1] Ye Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin 
      "PlatEMO: A MATLAB Platform for Evolutionary Multi-Objective Optimization [Educational Forum]" 
      IEEE Computational Intelligence Magazine, 2017, 12(4): 73-87.
  
  
  
  
  
  
