function [KeepV, KeepVR,KeepVP,KeepDR,KeepDP] = StoichiometryViaMILP2(Atomic,MaxNumReactants,MaxNumProducts,MaxMolReactants,MaxMolProducts,MaxLoops)
%
% Determine all feasible stoichiometies via MILP. Uses the Atomic matrix to
% define equality constaints, binary variables are used to define sucessive
% cuts, to remove previous solutions. The algorithm then iterates to find
% all feasible solutions.
%
% Written: Mark Willis (11/1/17)
% Last Update:    
%
% Chemical reactions are modelled as Vrx --> Vpx e.g (2,0,0,0)x --> (0 1 0
% 0)x gives V=Vp-Vr =(-2 1 0 0)
%

KeepV=[]; KeepVR=[];KeepVP=[]; KeepDR=[]; KeepDP=[];

[NumAtomicConstraints, NumSpecies]=size(Atomic);
%
% Define the total number of independent variables which are the stoichiometric
% numbers and associated binary variables

NumIndependentVariables = 4*NumSpecies;

%
% Set the MILP constraints
%

xvars= 1:2*NumSpecies;                                 % Stoichiometric numbers
xvarsBinary = 2*NumSpecies+1:NumIndependentVariables;  % Binary variables

% Set bounds on the model parameters
L=zeros(1,2*NumSpecies); U=2*ones(1,2*NumSpecies);
lb(xvars) = L;
ub(xvars) = U;    
lb(xvarsBinary) = zeros(2*NumSpecies,1);
ub(xvarsBinary) = ones(2*NumSpecies,1);

A=[]; b=[]; Aeq=[];Beq=[];

%Binary variable constraints
ABin1 = horzcat(eye(2*NumSpecies), -(zeros(2*NumSpecies)+diag(U')));
A = [A;ABin1];

%Binary variable constraints
ABin2 = horzcat(-eye(2*NumSpecies), (zeros(2*NumSpecies)+diag(L'+1)));
A = [A;ABin2];

% Fill up the 'b' vector where Ax <= b
b = [b;zeros(2*NumSpecies,1); zeros(2*NumSpecies,1)];

% Atomic matrix constraint, Atomic*x =0

Aeq = horzcat(-Atomic, Atomic, zeros(NumAtomicConstraints,2*NumSpecies));
Beq = zeros(NumAtomicConstraints,1);

% At least one species reacts and one is produced...

temp = horzcat(zeros(1,2*NumSpecies), -ones(1,NumSpecies),zeros(1,NumSpecies));
A=[A;temp];
b = [b;-1];

temp = horzcat(zeros(1,2*NumSpecies), zeros(1,NumSpecies),-ones(1,NumSpecies));
A=[A;temp];
b = [b;-1];

% The same species cant be a reactant and a product...

temp = horzcat(zeros(NumSpecies,2*NumSpecies), eye(NumSpecies),eye(NumSpecies));
A=[A;temp];
b = [b;ones(NumSpecies,1)];

% Specify maximum number of reactants and products

temp = horzcat(zeros(1,2*NumSpecies), zeros(1,NumSpecies),ones(1,NumSpecies));
A=[A;temp];
b = [b;MaxNumReactants];

temp = horzcat(zeros(1,2*NumSpecies), ones(1,NumSpecies),zeros(1,NumSpecies));
A=[A;temp];
b = [b;MaxNumProducts];

% Set maximum molecularity of reactants and products

temp = horzcat(zeros(1,NumSpecies),ones(1,NumSpecies),zeros(1,2*NumSpecies));
A=[A;temp];
b = [b;MaxMolReactants];

temp = horzcat(ones(1,NumSpecies),zeros(1,NumSpecies),zeros(1,2*NumSpecies));
A=[A;temp];
b = [b;MaxMolProducts];

% Define the objective function; minimise the sum of the stoichiometric
% numbers

f = [ones(2*NumSpecies,1); zeros(2*NumSpecies,1)];

% Call the MILP - note the additional constraint that the stoichiometic
% numbers are integers...

for i=1:MaxLoops    % Loop a set number of times...

options = optimoptions('intlinprog','Display','off');
[xLinInt,fval,exitFlagInt,output] = intlinprog(f,[xvars xvarsBinary],A,b,Aeq,Beq,lb,ub,options);

% Stop looping when no feasible solution is found
if exitFlagInt == -2
    break;
end

xpars =xLinInt(xvars);
xparsBinary= xLinInt(xvarsBinary);

% Stoichiometry is...

VR = xpars(1:NumSpecies);
VP = xpars(NumSpecies+1:2*NumSpecies);
V = VP-VR;

% Binary variables are....

DR = xparsBinary(1:NumSpecies);
DP = xparsBinary(NumSpecies+1:2*NumSpecies);

% Store all suggested stoichiometies

KeepV =[KeepV V];
KeepVR =[KeepVR VR];
KeepVP =[KeepVP VP];
KeepDR =[KeepDR DR];
KeepDP =[KeepDP DP];

% Define the integer cut...

IndexZeros=find(round(xparsBinary)==0);
IndexOnes=find(round(xparsBinary)==1);
NumberOfOnes=length(IndexOnes);
NumberOfZeros=length(IndexZeros);
temp=zeros(1,2*NumSpecies);
temp(IndexOnes)=1;
temp(IndexZeros)=-1;
temp2=2*NumSpecies-1-NumberOfZeros;

% Add the integer cut to the inequality constraints...

temp3 = horzcat(zeros(1,2*NumSpecies), temp);
A=[A;temp3];
b = [b;temp2];

end  % End MILP loop