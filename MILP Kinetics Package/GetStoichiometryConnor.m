%% Get stoichiometry.
% Uses reaction matrix to get all feasible stoichiometric terms
% The solution is via a integer linear program (ILP)

% Calls: StoichiometryViaMILP2
%
% Mark Willis & Connor Taylor

% clear all;
% close all;

%% Set the reaction example
%

ReactionExample = 3;

%% Define the algorithm control parameters

% No more than two reactant species and two product species

MaxNumReactants = 2; 
MaxNumProducts = 2; 

% Molecularity of reactants and products

MaxMolReactants = 2;
MaxMolProducts = 2;

%% Set the atomic matrix (example 1)
% For this example it is assumed there are 4 species and Atomic is a list
% of the relative molecular weights of the species

if ReactionExample == 1

    Atomic =[1 1 1 2];

end

%% Set the atomic matrix (example 1)
% For this example it is assumed there are 5 species. 
% Each column of Atomic is a species. Each row a chemical element (C, H, O)

if ReactionExample == 2
    
Atomic = [4 2 4 6 0;10 4 6 12 2;1 2 3 2 1];

end 

%% Atomic matrix CT added reaction, Multistep Conversion (C, H, N, O, Cl)
% A, B, C, D, E, F(Hcl)
if ReactionExample == 3
    
    Atomic = [4 4 8 8 12 0;2 9 10 10 18 1;2 1 3 3 4 0;0 1 1 1 2 0;2 0 1 1 0 1]
end


%% Set maximum number of loops and call the optimisation function
% KeepV is a list of all feasible stoichiometry (subject to constraints)

MaxLoops = 1000;

[KeepV, KeepVR,KeepVP,KeepDR,KeepDP] = StoichiometryViaMILP2(Atomic,MaxNumReactants,MaxNumProducts,MaxMolReactants,MaxMolProducts,MaxLoops);

KeepV
