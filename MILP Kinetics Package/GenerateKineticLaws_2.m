function [TotalCombinationsFunctions,NumberofKValues] = GenerateKineticLaws_2(KeepV,a,TotalCombinations,Combinations,Conc)
%Function used to assign all possible reaction models.

Counter = 0;
Alphabet = ['A';'B';'C';'D'];
OrderedRateLaws = [];
OrderedRateLaw = []; %used as a temporary variable.
TotalCombinationsFunctions = {'1','2','3','4','5','6','7','8','9','10'};
TotalCombinationsFunctions = TotalCombinationsFunctions(1:size(Conc,2));
NumberofKValues = [];

for i = 1:size(KeepV,2) %for each Combinations.CombX (X = 1,2,etc.).
    eval(strcat(sprintf('CombinationSize = size(Combinations.Comb%d,1);',i)))
    
    for L = 1:CombinationSize %for each combination of KeepV columns.
        clear KValues
        KValues = struct('Placeholder',1);
        eval(strcat(sprintf('CombColumns = Combinations.Comb%d(L,:);',i)))

NewKeepV = [];
for X = 1:size(CombColumns,2)
    NewKeepV(:,X) = KeepV(:,CombColumns(X)); %only the relevant reactions from KeepV.
end

EquationsNeg = [];
EquationsPos = [];
for M = 1:size(NewKeepV,2)
    EquationsNeg{M} = (sprintf('-KV.K%d',M)); %beginning of each equation for each column.
    EquationsPos{M} = (sprintf('KV.K%d',M));
end

Positives = NewKeepV;%matrix for concentrations going up.
Negatives = NewKeepV;%matrix for concentrations going down.


Positives(Positives < 0) = 0;
Negatives(Negatives > 0) = 0;
for M = 1:size(Negatives,2)%for each column.
    for N = 1:size(Negatives,1)%for each entry in the column.
        if Negatives(N,M) < 0 
            %assign a rate equation based on the negative values, meaning
            %that the rate laws are dictated by decreasing [A] or [B] etc.
            power = num2str(abs(Negatives(N,M)));
            EquationsNeg{M} = strcat(EquationsNeg{M},sprintf('*(a.%s^',Alphabet(N)),sprintf(power),sprintf(')'));
            EquationsPos{M} = strcat(EquationsPos{M},sprintf('*(a.%s^',Alphabet(N)),sprintf(power),sprintf(')'));
        else
        end
    end
end


for M = 1:size(NewKeepV,2)
    for N = 1:size(NewKeepV,1)
        if Positives(N,M) > 0
            %this makes sure that if there is a reaction such as A + B -> C
            %+ C, then d[C]/dt = 2*rate (to account for the double C
            %concentration).
            multiply = num2str(Positives(N,M));
            EquationsPos{M} = strcat(multiply,sprintf('*'),EquationsPos{M});
            break
        else
        end
    end
end


for N = 1:size(NewKeepV,1)
OrderedRateLaws{N} = '@(KV,a)0';
end

%now to look through each row to order the rate laws.
for N = 1:size(NewKeepV,1)
    for M = 1:size(NewKeepV,2)
        if NewKeepV(N,M) < 0
                OrderedRateLaws{N} = strcat(OrderedRateLaws{N},sprintf('+'),EquationsNeg{M});
        elseif NewKeepV(N,M) > 0
                OrderedRateLaws{N} = strcat(OrderedRateLaws{N},sprintf('+'),EquationsPos{M});
        else
    end
    end
end


TotalCombinations;
Functions = [];
for functioncounter = 1:size(OrderedRateLaws,2)
    Functions{functioncounter} = str2func([OrderedRateLaws{functioncounter}]);
end

Counter = Counter+1;
 TotalCombinationsFunctions(Counter,:) = Functions;
NumberofKValues(Counter) = size(CombColumns,2);
sprintf('Number of functions added: %d of %d',Counter,TotalCombinations)
    end
end
Counter


