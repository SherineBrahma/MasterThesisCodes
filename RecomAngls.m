function RecomAnglArry = RecomAngls(G, NoOfAnglUsed, AnglArry, NoOfTimeInstances )
%{
clear all
NoOfAnglUsed = 3;
AnglArry = [10 20 30 40];
NoOfTimeInstances = 3;
%}

TotlAngl = size(AnglArry,2);

disp('....Determining Threshold....')
%{
SampG = randsample(reshape(G,[],1),1000000) ;
%hist(SampG,10000)
[Mean, Std] = GaussFit(SampG);
Mean
Std
clear SampG
Trhld = Mean + 3*Std;              %User Specified thersholding for filtering out smaller angles
%}
Trhld = 1; % No Thresholding

disp('.... Computing vector independency indication matrix....')

for URowIndx = 1:1:TotlAngl
    for UColIndx = 1:1:TotlAngl
        SubGThdRow = (URowIndx-1) * NoOfTimeInstances + 1:1:(URowIndx-1) * NoOfTimeInstances + 1 + (NoOfTimeInstances-1);        
        SubGThdCol = (UColIndx-1) * NoOfTimeInstances + 1:1:(UColIndx-1) * NoOfTimeInstances + 1 + (NoOfTimeInstances-1);        
        SubG = G(SubGThdRow, SubGThdCol);
        SubGThd = SubG;
        SubGThd(SubG > Trhld) = 1;
        U(URowIndx, UColIndx) = sum(sum(SubGThd));
    end
end

% Make the diagonals 0
U = U - diag(diag(U));

disp('.... Getting recommended angles....')

umin = sum(sum(U));
SrcIndx = 1:1:size(AnglArry,2);
uminarry = [];
for ItCnt = 1:1:1000000
    RndPerm = randperm(size(AnglArry,2));
    x = sort(SrcIndx(RndPerm(1:(NoOfAnglUsed))));
    u = sum(sum(U(x,x)));
    uminarry(ItCnt) = umin;
    if u < umin
        umin = u;
        xmin = x;
    end
end

plot(uminarry(2:1:end))
RecomAnglArry = AnglArry(xmin);

end