function RecomAnglArry = RecomAngls(SysMat, NoOfAnglUsed, AnglArry, NoOfTimeInstances )
%{
clear all
SysMat = [ 1 2 3   4  5  6  1 2 3    7 6 2;
           7 8 9  11 27 23  7 8 9    2 4 1;
           3 4 1  34 22 11  3 4 1   13 17 15];
SysMat = [ 1 2 3   4  5  6  1 2 3;
           7 8 9  11 27 23  7 8 9;
           3 4 1  34 22 11  3 4 1];       
NoOfAnglUsed = 3;
AnglArry = [10 20 30 40];
NoOfTimeInstances = 3;
%}
TotlAngl = size(AnglArry,2);
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
G = abs(G);
%B = single(G);
%gpuG = gpuArray(B);
%a = reshape(G,[],1);

y = randsample(a,100000) ;

[Mean, Std] = GaussFit(reshape(B,[],1));

plot(reshape(G,[],1))

disp('.... Computing vector independency indication matrix....')

for URowIndx = 1:1:TotlAngl
    for UColIndx = 1:1:TotlAngl
        display(URowIndx)
        display(UColIndx)
        SubGThdRow = (URowIndx-1) * NoOfTimeInstances + 1:1:(URowIndx-1) * NoOfTimeInstances + 1 + (NoOfTimeInstances-1);        
        SubGThdCol = (UColIndx-1) * NoOfTimeInstances + 1:1:(UColIndx-1) * NoOfTimeInstances + 1 + (NoOfTimeInstances-1);        
        SubG = abs(G(SubGThdRow, SubGThdCol));
        SubGThd = SubG;
        %gpuSubGThd = gpuArray(reshape(SubG,[],1));
        [Mean, Std] = GaussFit(reshape(SubG,[],1));
        Trhld = Mean + 4*Std; % mean(mean(SubG)); %User Specified thersholding for filtering out smaller angles
        SubGThd(SubG < Trhld) = 1;
        SubGThd(SubG >= Trhld) = 0;
        U(URowIndx, UColIndx) = sum(sum(SubGThd));
    end
end

%hist(reshape(SubG,[],1),516)
% Make the diagonals 0
U = U - diag(diag(U));

disp('.... Getting recommended angles....')
%NoOfAnglUsed = 36;
% To get the best angles
[Usor, UsorIndx] = sort(U,2,'descend');
Ucum = cumsum(Usor,2);                                                % Accumulate row values
[MaxCumIP, MaxCumIPIndx] = max(Ucum(:,NoOfAnglUsed));                 % Choose the best value for corresponding number of angles
if NoOfAnglUsed <= 1
    disp('error - No of Angle used should be atleast 2');
    beep;    return
end
AIndx = [MaxCumIPIndx UsorIndx(MaxCumIPIndx, 1:(NoOfAnglUsed - 1))];  % Get the angles which are used for the best angles     
RecomAnglArry = sort(AnglArry(AIndx));

end











