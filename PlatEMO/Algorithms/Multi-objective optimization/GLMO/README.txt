
 Copyright (C) 2020 Heiner Zille

 This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 
 International License. (CC BY-NC-SA 4.0). To view a copy of this license, 
 visit http://creativecommons.org/licenses/by-nc-sa/4.0/ or see the 
 pdf-file "License-CC-BY-NC-SA-4.0.pdf" that came with this code. 

 You are free to: 
 * Share ? copy and redistribute the material in any medium or format
 * Adapt ? remix, transform, and build upon the material 
 Under the following terms:
 * Attribution ? You must give appropriate credit, provide a link to the 
    license, and indicate if changes were made. You may do so in any reasonable 
    manner, but not in any way that suggests the licensor endorses you or your use.
 * NonCommercial ? You may not use the material for commercial purposes.
 * ShareAlike ? If you remix, transform, or build upon the material, you must 
   distribute your contributions under the same license as the original.
 * No additional restrictions ? You may not apply legal terms or technological 
   measures that legally restrict others from doing anything the license permits.

 Author of this Code: 
  Heiner Zille <heiner.zille@ovgu.de> or <heiner.zille@gmail.com>

 This code is based on the following publications:

 1) Heiner Zille 
    "Large-scale Multi-objective Optimisation: New Approaches and a Classification of the State-of-the-Art"  
    PhD Thesis, Otto von Guericke University Magdeburg, 2019 
    http://dx.doi.org/10.25673/32063 

 2) Heiner Zille, Hisao Ishibuchi, Sanaz Mostaghim and Yusuke Nojima
    "Mutation Operators Based on Variable Grouping for Multi-objective Large-scale Optimization"
    IEEE Symposium Series on Computational Intelligence (SSCI), IEEE, Athens, Greece, December 2016
    https://ieeexplore.ieee.org/document/7850214 

 This file is intended to work with the PlatEMO framework version 2.5. 
 Date of publication of this code: 06.04.2020 
 Last Update of this code: 06.04.2020 
 A newer version of this algorithm may be available. Please contact the author 
 or see http://www.ci.ovgu.de/Research/Codes.html. 
  
 To use this program, simply copy the complete "GLMO"-folder and its contained files 
 into the "Algorithms"-folder of the PlatEMO framework. GLMO should then be usable as any other 
 algorithm inside the framework. 
  
  File list:       
  	- GLMO/README.txt
	- GLMO/License-CC-BY-NC-SA-4.0.pdf
	
  	- GLMO/GLMO.m
  	- GLMO/GLMO_createGroups.m
  	- GLMO/GLMO_GA.m
  	- GLMO/GLMO_NSGAII.m
  	- GLMO/GLMO_NSGAIICrowdingDistance.m
  	- GLMO/GLMO_NSGAIIEnvironmentalSelection.m
  	- GLMO/GLMO_NSGAIII.m
  	- GLMO/GLMO_NSGAIIIEnvironmentalSelection.m
  	- GLMO/GLMO_SMPSO.m
  	- GLMO/GLMO_SMPSOOperator.m
  	- GLMO/GLMO_UpdateGbest.m
  	- GLMO/GLMO_UpdatePbest.m

The files may have been modified in Feb 2021 by the authors of the Platemo framework to work with the Platemo 3.0 release.   
