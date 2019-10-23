function [Error] = Reaction_2(input,KV, Results2, KeepV, Conc,a,CurrentCombinationsFunctions,I,TotalCombinations)
%Function to be optimised that calls ODE solver and outputs 'Error' as the
%variable to be optimised. K values are iterated to find the minimum sum of
%squared error between the simulated data and the experimental data.

%Assigns all of the inputs from the optimiser.
for X = 1:size(input,2)
    eval(sprintf('KV.K%d = input(X);',X))
end

%Times in seconds for the ODE solver to evaluate.
Time = [0,15,30,45,60,120,240,360,480,600];

interupt_time = 400; %Maximum time in seconds for the ODE solver to run.

%Sets the output function to incorporate interuption time:
outputFun = @(Time,Conc,flag)interuptFunODE(Time,Conc,flag,interupt_time);
opts = odeset('OutputFcn',outputFun);
warning off

%This is the maximum allowed total molar value of all species, this is
%defined so that ODEs don't form a 'runaway integration' where values
%quickly escalate by orders of magnitude.
TotalConc = sum(Conc)*10;


    %ODE solver:
    [TimeData,ConcData] = ode15s(@(Time,Conc)ReactionKineticLaws_2(Time,Conc,KV,a,CurrentCombinationsFunctions,I,TotalConc),Time,Conc,opts);
 
    %Sum of squared error comparison, serves as a comparison to
    %experimental results after the particular run using the particular K
    %inputs.
    Error = sse(ConcData-Results2);  




