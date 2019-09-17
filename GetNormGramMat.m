function G = GetNormGramMat(SysMat)
%{
clear all
SysMat = [ 1 2 3   4  5  6  1 2 3    7 6 2;
           7 8 9  11 27 23  7 8 9    2 4 1;
           3 4 1  34 22 11  3 4 1   13 17 15];
SysMat = [ 1 2 3   4  5  6  1 2 3;
           7 8 9  11 27 23  7 8 9;
           3 4 1  34 22 11  3 4 1];
%}
A = SysMat'; 
clear SysMat

% To normalize the columns
nn = sqrt(sum(A.*conj(A),1));
if ~all(nn)
    disp('error - input contains a zero column');
    u=NaN;   beep;    return
end
nA = bsxfun(@rdivide,A,nn);  % nA is a matrix with normalized columns

disp('.... Calculating Gram matrix....')

% Get the Gram matrix

G = (nA')*nA;
clear A nA
G = abs(G);

%{
ChnkFact = 8;
ChnkSize = ChnkFact*NoOfTimeInstances;  % Computing Absolute in chunks to save memory
AbsG = [];
for AbsGRowIndx = 1:1:TotlAngl/ChnkFact
    AbsSubGCol = [];
    display(AbsGRowIndx)
    for AbsGColIndx = 1:1:TotlAngl/ChnkFact
        SubGRow = (AbsGRowIndx-1) * ChnkSize + 1:1:(AbsGRowIndx-1) * ChnkSize + 1 + (ChnkSize-1);        
        SubGCol = (AbsGColIndx-1) * ChnkSize + 1:1:(AbsGColIndx-1) * ChnkSize + 1 + (ChnkSize-1);
        IntAbsSubG = abs(G(SubGRow, SubGCol));
        AbsSubGCol = [AbsSubGCol IntAbsSubG];
    end
    AbsG = [AbsG; AbsSubGCol];
end
G = AbsG;
clear AbsG;
%}

end











