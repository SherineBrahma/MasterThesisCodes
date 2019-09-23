function RecomAnglArry = RecomAnglsHor(SysCell, NoOfAnglUsed, AnglArry, NoOfTimeInstances )
%{
clear all
NoOfAnglUsed = 3;
AnglArry = [10 20 30 40];
NoOfTimeInstances = 3;
%}

TotlAngl = size(AnglArry,2);

 SysMat = horzcat(SysCell{:,1});
 
for URowIndx = 1:1:TotlAngl
    SubS1Row = 1:1:NoOfTimeInstances;        
    SubS1Col = (URowIndx-1) * 4225 + 1:1:(URowIndx-1) * 4225 + 1 + (4225-1);        
    SubS1 = SysMat(SubS1Row, SubS1Col);
    gpuSubS1 = gpuArray(SubS1);
    nn1 = sqrt(sum(gpuSubS1.*conj(gpuSubS1),1));
    NorS1 = bsxfun(@rdivide,gpuSubS1,nn1);
    %GSubPivot = abs((NorS1')*NorS1);
    for UColIndx = 1:1:TotlAngl
        SubS2Row = 1:1:NoOfTimeInstances;        
        SubS2Col = (UColIndx-1) * 4225 + 1:1:(UColIndx-1) * 4225 + 1 + (4225-1);        
        SubS2 = SysMat(SubS2Row, SubS2Col);
        gpuSubS2 = gpuArray(SubS2);
        nn2 = sqrt(sum(gpuSubS2.*conj(gpuSubS2),1));
        NorS2 = bsxfun(@rdivide,gpuSubS2,nn2);
        %GSubLoop = abs((NorS2')*NorS2);
        GSubLoop = abs((NorS1')*NorS2);
        U(URowIndx, UColIndx) = gather(det(GSubLoop));
    end
end

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


