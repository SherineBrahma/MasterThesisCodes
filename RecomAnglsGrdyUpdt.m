function [TimeElpsd, RecomAnglArry, cnminArry] = RecomAnglsGrdyUpdt(SysCell, NoOfAnglUsed, AnglArry)
    tic
    SysCell(:,2) = cellfun(@(sysmat)  cond(sysmat), SysCell, 'UniformOutput',false);

    NoOfAngl = size(AnglArry,2);
    [~,InitIndx] =  min(cell2mat(SysCell(:,2)));
    SysMat = SysCell{InitIndx,1};
    if size(SysMat,2) > size(SysMat,1)
        PadSize = size(SysMat,2)-size(SysMat,1);
        SysMat = padarray(SysMat,[size(SysMat,2)-size(SysMat,1) 0],0,'pre');
        PrePadFlag = 1;
    else
        PrePadFlag = 0;
    end
    [U,S,V]=svd(SysMat,0);

    AnglIndx = [InitIndx];

    SrcIndx = 1:1:NoOfAngl;
    SrcIndx(SrcIndx==InitIndx)=[];
    cnmin = inf;
    ULOutLoop = U;

disp('....Employing greedy approach to find angles....')
    for AnglUsedIndx = 1:1:(NoOfAnglUsed-1)
        for MatIndx = 1:1:size(SrcIndx,2)
            LoopMat = SysCell{SrcIndx(MatIndx),1};
            UInLoop = [ULOutLoop; zeros(size(LoopMat))];
            A = [zeros(size(SysMat,1),size(LoopMat,1)); eye(size(LoopMat,1))]; 
            B = LoopMat';
            fprintf('Angles remaining to be searched: %d\n', size(SrcIndx,2)-MatIndx + 1);
            [UInLoop,SInLoop,V,time] = updatesvd(UInLoop,S,V,A,B);
            fprintf('Time per search: %d\n', time);
            dS = diag(SInLoop); 
            cn(MatIndx) = dS(1)/dS(end);
            if cn(MatIndx) < cnmin
               cnmin = cn(MatIndx);
               cnminArry(MatIndx) = cnmin;
               MinAnglIndx = SrcIndx(MatIndx);
            else
                cnminArry(MatIndx) = cnmin;
            end
        end
        cnmin = inf;
        ULOutLoop = UInLoop;
        SysMat = [SysMat; SysCell{MinAnglIndx,1}];
        AnglIndx = [AnglIndx MinAnglIndx];
        SrcIndx(SrcIndx==MinAnglIndx)=[];
        fprintf('Angles found so far: %d\n', AnglUsedIndx);
    end

    if PrePadFlag == 1
        SysMat(1:1:PadSize,:) = [];
        PrePadFlag = 1;
    end

    RecomAnglArry = AnglArry(sort(AnglIndx));
    TimeElpsd = toc;
end
