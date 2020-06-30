function Population = EnvironmentalSelection(Population,Offspring,Z,nDivs)
% The environmental selection of MOEA/PC

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Roman Denysiuk

% for offspring, determine grid index and polar coordinates
[off.idx, off.theta, off.rho] = PolarCoordinatesGrid(Offspring.obj, Z, nDivs);

% for population member corresponding to grid index of offspring,
% determine true grid index and polar coordinates
[pop.idx, pop.theta, pop.rho] = PolarCoordinatesGrid(Population(off.idx).obj, Z, nDivs);

% update grid
if (off.idx == pop.idx) 
    % pop member and off reside in the same grid     
    % offspring replaces pop member if it has a smaller radius
    if off.rho < pop.rho
        Population(off.idx) = Offspring; 
    end    
    return
    
elseif any(Offspring.obj < Population(off.idx).obj)
    % pop memebr and offspring are swapped because
    % pop member resides in other grid 
    tmp = Population(off.idx);
    Population(off.idx) = Offspring;
    Offspring = tmp;
    Population = EnvironmentalSelection(Population,Offspring,Z,nDivs);
    
end

end