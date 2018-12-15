function value = WOF_transformationFunction(xPrime,weight,maxVal,minVal,method)
% Implements the transformation functions used in WOF-SMPSO. Three 
% methods can be chosen. The p-Value and Multiplication-transformations
% have been introduced in publication (2), see above. 
% The Intervall-transformation (parameter free) has been introduced in 
% publication (2), see above. The interval-intersection method from (2)
% is currently not implemented. 

% ----------------------------------------------------------------------- 
%  WOF_transformationFunction.m 
%  Copyright (C) 2018 Heiner Zille
% 
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%  Author of this Code: 
%   Heiner Zille <heiner.zille@ovgu.de>
%
%  This file belongs to the following publications:
%
%  1) Heiner Zille and Sanaz Mostaghim
%     "Comparison Study of Large-scale Optimisation Techniques on the LSMOP Benchmark Functions"  
%     IEEE Symposium Series on Computational Intelligence (SSCI), IEEE, Honolulu, Hawaii, November 2017
%     https://ieeexplore.ieee.org/document/8280974 
% 
%  2) Heiner Zille, Hisao Ishibuchi, Sanaz Mostaghim and Yusuke Nojima
%     "A Framework for Large-scale Multi-objective Optimization based on Problem Transformation"
%     IEEE Transactions on Evolutionary Computation, Vol. 22, Issue 2, pp. 260-275, April 2018.
%     http://ieeexplore.ieee.org/document/7929324
%  
%  3) Heiner Zille, Hisao Ishibuchi, Sanaz Mostaghim and Yusuke Nojima
%     "Weighted Optimization Framework for Large-scale Mullti-objective Optimization"
%     Genetic and Evolutionary Computation Conference (GECCO), ACM, Denver, USA, July 2016
%     http://dl.acm.org/citation.cfm?id=2908979
%
%  Date of publication: 12.10.2018 
%  Last Update: 12.10.2018
% -----------------------------------------------------------------------

    value = xPrime;
    switch method
        case 1 %multiplication
            value = xPrime*weight;
        case 2 %p-value
            pWert = 0.2;
            value = xPrime+pWert*(weight-1.0)*(maxVal-minVal);
        case 3 %interval
            if weight > 1.0
                weight = weight - 1.0;
                interval = maxVal - xPrime;
                value = xPrime + weight * interval;
            else
                interval = xPrime - minVal;
                value = minVal + weight * interval;
            end           
    end
    
    %do repair
    if value < minVal
       value = minVal;
    elseif value > maxVal
       value = maxVal;
    end
    
end