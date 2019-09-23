clear all
close all
 
SysMat1 = [11 56 ;
           3 28;
           78 31;
           56 25];
SysMat2 = [ 42 -9 ;
           78 31;
           25 75;
           53 89];
SysMat3 = [ 64 13;
           30 34;
           30 34;
           42 11];
SysMat4 = [ 31 10;
           30 34;
           23 43;
           24 19];
SysMat5 = [ 12 74;
           23 43;
           64 13;
           19 12];
       
SysCell = {SysMat1; SysMat2; SysMat3; SysMat4; SysMat5};
SysCell(:,2) = cellfun(@(sysmat)  cond(sysmat), SysCell, 'UniformOutput',false);

NoOfAngl = 5;
NoOfAnglUsd = 2;
NoOfPixels = 2;
[~,InitIndx] =  min(cell2mat(SysCell(:,2)));
SysMat = SysCell{InitIndx,1};
[U,S,V]=svd(SysMat,0);

AnglIndx = [InitIndx];

SrcIndx = 1:1:NoOfAngl;
SrcIndx(SrcIndx==InitIndx)=[];
cnmin = inf;
ULOutLoop = U;

for AnglUsdIndx = 1:1:(NoOfAnglUsd-1)
    for MatIndx = 1:1:size(SrcIndx,2)
        LoopMat = SysCell{SrcIndx(MatIndx),1};
        UInLoop = [ULOutLoop; zeros(size(LoopMat))];
        A = [zeros(size(SysMat,1),size(LoopMat,1)); eye(size(LoopMat,1))]; 
        B = LoopMat';
        [UInLoop,SInLoop,V,time] = updatesvd(UInLoop,S,V,A,B);
        dS = diag(SInLoop); 
        cn(MatIndx) = dS(1)/dS(end);
        if cn(MatIndx) < cnmin
           cnmin = cn(MatIndx)
           MinAnglIndx = SrcIndx(MatIndx);
        end
    end
    cnmin = inf;
    ULOutLoop = UInLoop;
    SysMat = [SysMat; SysCell{SrcIndx(MinAnglIndx),1}];
    AnglIndx = [AnglIndx MinAnglIndx];
    SrcIndx(SrcIndx==MinAnglIndx)=[];
end

%{
ImpSysMat = [SysMat1; SysMat2]; 
c12 = cond(ImpSysMat)
ImpSysMat = [SysMat2; SysMat5]; 
c23 = cond(ImpSysMat)
ImpSysMat = [SysMat1; SysMat3; SysMat5]; 
c13 = cond(ImpSysMat)
%}
