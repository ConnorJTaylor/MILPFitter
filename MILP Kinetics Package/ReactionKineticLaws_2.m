function [dydt] = ReactionKineticLaws_2(Time,Conc,KV,a,CurrentCombinationsFunctions,I,TotalConc)
%This function is the ODE function where the K values are unpacked and the
%rate laws are used to see how concentrations change over time.

dydt = zeros(size(Conc));

%Assigns alphabet terms to the each of the concentrations.
try 
    A = Conc(1);
    a.A = A;
    B = Conc(2);
    a.B = B;
    C = Conc(3);
    a.C = C;
    D = Conc(4);
    a.D = D;
    E = Conc(5);
    a.E = E;
    F = Conc(6);
    a.F = F;
catch
end

%Generate an error if the ODE solver begins a 'runaway integration',
%therefore stopping the integration:
% if sum(Conc) > TotalConc || isnan(sum(Conc)) 
%     error
% end

%Evaluate dydt:
for M = 1:size(Conc,1) %for each model term.
    %Assign fHandle as the current reaction model term:
    fHandle = CurrentCombinationsFunctions{I,M};
    %Evaluate current reaction model term:
    dydt(M) = fHandle(KV,a);
end