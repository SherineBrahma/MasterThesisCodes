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

disp('.... Computing vector independency indication matrix....')

for URowIndx = 1:1:TotlAngl
    for UColIndx = 1:1:TotlAngl
        SubGRow = (URowIndx-1) * NoOfTimeInstances + 1:1:(URowIndx-1) * NoOfTimeInstances + 1 + (NoOfTimeInstances-1);        
        SubGCol = (UColIndx-1) * NoOfTimeInstances + 1:1:(UColIndx-1) * NoOfTimeInstances + 1 + (NoOfTimeInstances-1);        
        SubG = G(SubGRow, SubGCol);
        gpuSubG = gpuArray(SubG);
        SorSubG = sort(reshape(gpuSubG,[],1),'ascend');
        NoColConsd = round(0.6*size(SorSubG,1));
        U(URowIndx, UColIndx) = gather(sum(sum(SorSubG(1:1:NoColConsd))));
    end
end
%load('UMinEig','U');

%save('CohMat72Ang','U')
%save('UMinEig','U');
disp('.... Getting recommended angles....')
% 0-360: 7.1668, 83.9017   My: 4.1420, 44.7674
NoOfAnglUsed = 9;
umin = sum(sum(U));
SrcIndx = 1:1:size(AnglArry,2);
uminarry = [];
for ItCnt = 1:1:1000000
    RndPerm = randperm(72);%size(AnglArry,2));
    x = sort(SrcIndx(RndPerm(1:(NoOfAnglUsed))));
    %x = sort(SrcIndx(1:2:72));
    %x = [intx (intx + 18*ones(size(intx))) (intx + 2*18*ones(size(intx))) (intx + 3*18*ones(size(intx)))];
    u = sum(sum(U(x,x)));
    uminarry(ItCnt) = umin;
    if u <= umin
        umin = u;
        xmin = x;
    end
end

%plot(uminarry(2:1:end))
RecomAnglArry = AnglArry(xmin);

end


