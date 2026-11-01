function Probabilities = FLCS(RewardGA, RewardDE, Diversity, Convergence)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Yubo Wang (email: 1202111230@cug.edu.cn)

    persistent FuzzySystem
    if isempty(FuzzySystem)
        FuzzySystem = mamfis('Name', 'Operator Selection');
        FuzzySystem = addInput(FuzzySystem, [0 1], 'Name', 'RewardGA');
        FuzzySystem = addInput(FuzzySystem, [0 1], 'Name', 'RewardDE');
        FuzzySystem = addInput(FuzzySystem, [0 1], 'Name', 'Diversity');
        FuzzySystem = addInput(FuzzySystem, [0 1], 'Name', 'Convergence');
        FuzzySystem = addOutput(FuzzySystem, [0 1], 'Name', 'GA_DoM');
        FuzzySystem = addOutput(FuzzySystem, [0 1], 'Name', 'DE_DoM');
        FuzzySystem = addMF(FuzzySystem, 'RewardGA', 'trapmf', [0 0 0.3 0.4], 'Name', 'Low');
        FuzzySystem = addMF(FuzzySystem, 'RewardGA', 'trimf', [0.3 0.5 0.7], 'Name', 'Mid');
        FuzzySystem = addMF(FuzzySystem, 'RewardGA', 'trapmf', [0.6 0.7 1 1], 'Name', 'High');
        FuzzySystem = addMF(FuzzySystem, 'RewardDE', 'trapmf', [0 0 0.3 0.4], 'Name', 'Low');
        FuzzySystem = addMF(FuzzySystem, 'RewardDE', 'trimf', [0.3 0.5 0.7], 'Name', 'Mid');
        FuzzySystem = addMF(FuzzySystem, 'RewardDE', 'trapmf', [0.6 0.7 1 1], 'Name', 'High');
        FuzzySystem = addMF(FuzzySystem, 'Diversity', 'trapmf', [0 0 0.3 0.4], 'Name', 'Low');
        FuzzySystem = addMF(FuzzySystem, 'Diversity', 'trimf', [0.3 0.5 0.7], 'Name', 'Mid');
        FuzzySystem = addMF(FuzzySystem, 'Diversity', 'trapmf', [0.6 0.7 1 1], 'Name', 'High');
        FuzzySystem = addMF(FuzzySystem, 'Convergence', 'trapmf', [0 0 0.3 0.4], 'Name', 'Low');
        FuzzySystem = addMF(FuzzySystem, 'Convergence', 'trimf', [0.3 0.5 0.7], 'Name', 'Mid');
        FuzzySystem = addMF(FuzzySystem, 'Convergence', 'trapmf', [0.6 0.7 1 1], 'Name', 'High');
        FuzzySystem = addMF(FuzzySystem, 'GA_DoM', 'trapmf', [0 0 0.1 0.3], 'Name', 'Low');
        FuzzySystem = addMF(FuzzySystem, 'GA_DoM', 'trimf', [0.2 0.5 0.7], 'Name', 'Mid');
        FuzzySystem = addMF(FuzzySystem, 'GA_DoM', 'trapmf', [0.7 0.9 1 1], 'Name', 'High');
        FuzzySystem = addMF(FuzzySystem, 'DE_DoM', 'trapmf', [0 0 0.1 0.3], 'Name', 'Low');
        FuzzySystem = addMF(FuzzySystem, 'DE_DoM', 'trimf', [0.2 0.5 0.7], 'Name', 'Mid');
        FuzzySystem = addMF(FuzzySystem, 'DE_DoM', 'trapmf', [0.7 0.9 1 1], 'Name', 'High');
        ruleList = ["RewardGA = High & Convergence = Low  =>  GA_DoM = High"
                    "RewardGA = High & Convergence = Mid  =>  GA_DoM = Mid"
                    "RewardGA = High & Convergence = High  =>  GA_DoM = Mid"
                    "RewardGA = Mid & Convergence = Low  =>  GA_DoM = High"
                    "RewardGA = Mid & Convergence = Mid  =>  GA_DoM = Mid"
                    "RewardGA = Mid & Convergence = High  =>  GA_DoM = Low"
                    "RewardGA = Low & Convergence = Low  =>  GA_DoM = Mid"
                    "RewardGA = Low & Convergence = Mid  =>  GA_DoM = Low"
                    "RewardGA = Low & Convergence = High  =>  GA_DoM = Low"
                    "RewardDE = High & Diversity = Low => DE_DoM = High"
                    "RewardDE = High & Diversity = Mid => DE_DoM = Mid"
                    "RewardDE = High & Diversity = High => DE_DoM = Mid"
                    "RewardDE = Mid & Diversity = Low => DE_DoM = High"
                    "RewardDE = Mid & Diversity = Mid => DE_DoM = Mid"
                    "RewardDE = Mid & Diversity = High => DE_DoM = Low"
                    "RewardDE = Low & Diversity = Low => DE_DoM = Mid"
                    "RewardDE = Low & Diversity = Mid => DE_DoM = Low"
                    "RewardDE = Low & Diversity = High => DE_DoM = Low"];
        FuzzySystem = addRule(FuzzySystem, ruleList);
        FuzzySystem.DefuzzMethod = 'centroid';
    end
    inputs        = [RewardGA RewardDE Diversity Convergence];
    Probabilities = evalfis(FuzzySystem, inputs);
end