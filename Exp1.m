SysMat1 = round(rand([2,3])*100); %[11 56 43;
           %3 28 36];
SysMat2 = round(rand([2,3])*100);%[ 42 -9 56;
           %78 31 3];
SysMat3 = round(rand([2,3])*100);%[ 31 73 28;
          % 39 68 52];
SysMat = [SysMat1 SysMat2 SysMat3];
nn = sqrt(sum(SysMat.*conj(SysMat),1));
nA = bsxfun(@rdivide, SysMat, nn);
G = (nA')*nA;
Pixels = 3;
Angles = 3;
NoOfTimeInstances = 2;
for URowIndx = 1:1:Angles
    SubS1Row = 1:1:NoOfTimeInstances;        
    SubS1Col = (URowIndx-1) * Pixels + 1:1:(URowIndx-1) * Pixels + 1 + (Pixels-1);        
    SubS1 = SysMat(SubS1Row, SubS1Col);
    nn1 = sqrt(sum(SubS1.*conj(SubS1),1));
    NorS1 = bsxfun(@rdivide,SubS1,nn1);
    GSubPivot = abs((NorS1')*NorS1);
    for UColIndx = 1:1:Angles
        SubS2Row = 1:1:NoOfTimeInstances;        
        SubS2Col = (UColIndx-1) * Pixels + 1:1:(UColIndx-1) * Pixels + 1 + (Pixels-1);        
        SubS2 = SysMat(SubS2Row, SubS2Col);
        nn2 = sqrt(sum(SubS2.*conj(SubS2),1));
        NorS2 = bsxfun(@rdivide,SubS2,nn2);
        GSubLoop = abs((NorS2')*NorS2);
        GSub = abs((SubS1')*SubS2);
        %GsubNorm = abs((NorS1')*NorS2);
        GSubInvs(SubS1Col, SubS2Col) = (GSubPivot-GSubLoop);
        SorGSubInvs = sort(reshape(GSubInvs(SubS1Col, SubS2Col),1,[]),'descend');
        NoOfEntyConsd = round(1 * size(SorGSubInvs,2));
        U(URowIndx, UColIndx) = gather(sum(sum(SorGSubInvs(1:1:NoOfEntyConsd))));
        V(URowIndx, UColIndx) = gather(norm(SubS1-SubS2));
        %X(URowIndx, UColIndx) = gather(sum(diag(GsubNorm)));
    end
end

u(1) = -U(1,2) + U(2,1);
u(2) = -U(2,3) + U(3,2);
u(3) = -U(1,3) + U(3,1);

ImpSysMat = [SysMat1; SysMat2]; 
c12 = cond(ImpSysMat)
ImpSysMat = [SysMat2; SysMat3]; 
c23 = cond(ImpSysMat)
ImpSysMat = [SysMat1; SysMat3]; 
c13 = cond(ImpSysMat)

u/norm(u)
num2str(V)


