function [ReconsImgHTV, infoHTV, ReconsImgHBLSQR, infoHBLSQR]  = IRTikhRecon(SysMat, F, DataVec, OrigVec)
    
    %% Reshape Original Image
    OrigVec = reshape(OrigVec,[],1);

    %% Image Dimension
    ImgDim = [65,65];
    %acosd(.7265)
    
    %% Experiment Discrep of ICGLS &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    %{
    LambdaArry = [0.001];
    x = norm(DataVec-(SysMat*OrigVec))/ norm(DataVec);
    y = norm(SysMat'*(DataVec-(SysMat*OrigVec)))/ norm(SysMat'*DataVec);
    NoiseLevelCGLS = [x];
    for i = 1:1:size(NoiseLevelCGLS,2)
    optionsCGLS.RegParam = 'discrep';
    optionsCGLS.x_true = OrigVec;
    optionsCGLS.NoiseLevel = x;
    optionsCGLS.eta = 1.01;
    optionsCGLS.RegMatrix = F;
    optionsCGLS.NE_Rtol = 10^(-17);
    optionsCGLS.MaxIter = 500;
    %dbstop in IRcgls
    [TikhImageHBLSQR, infoCGLS] = IRcgls(SysMat, DataVec, [1:500], optionsCGLS);
    TikhImageCGLSStore = TikhImageHBLSQR;
    for j = 1:1:size(TikhImageHBLSQR,2)
    TikhImageHBLSQR = abs(reshape(TikhImageCGLSStore(:,j), ImgDim));
    ReconsImg = abs(TikhImageHBLSQR);
    ReconsImg = abs(reshape(ReconsImg, ImgDim));
    imagesc(ReconsImg)
    pause(.01);
    end
    end
 %}   
 %% Experiment Discrep of IRhtv &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
 
 x = norm(DataVec-(SysMat*OrigVec))/ norm(DataVec);
 y = norm(SysMat'*(DataVec-(SysMat*OrigVec)))/ norm(SysMat'*DataVec);
 NoiseLevelHTV = [x];
 %for i = 1:1:size(NoiseLevelHTV,2)
 optionsHTV.RegParam = 'discrep';%0;%'WGCV';%'WGCV';%
 %optionsHTV.GCVweight = 0;%'adapt';%
 optionsHTV.x_true = OrigVec;
 optionsHTV.inSolver = 'gmres';
 optionsHTV.adaptConstr = 'tvnn';
 optionsHTV.nonnegativity = 'on';
 optionsHTV.MaxIterIn = 30;
 optionsHTV.MaxIterOut = 30;
 optionsHTV.NoiseLevel = NoiseLevelHTV(1);%0;
 optionsHTV.eta = 1.01;
 %dbstop in IRhtv
 tic
 [TikhImageHTV, infoHTV] = IRhtv(SysMat, DataVec, optionsHTV);
 toc
 infoHTV.TimeElaspsed = toc;
 ReconsImgHTV = abs(TikhImageHTV);
 ReconsImgHTV = abs(reshape(ReconsImgHTV, ImgDim));
 %figure, imagesc(ReconsImgHTV)
 %end
 %options  = IRhtv('defaults');
 %{
    GCVNoiseEst = norm(DataVec-(SysMat*TikhImageHTV))/ norm(DataVec);
 %}
  %% Experiment Discrep of IRhybrid_lsqr &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  
    x = norm(DataVec-(SysMat*OrigVec))/ norm(DataVec);
    y = norm(SysMat'*(DataVec-(SysMat*OrigVec)))/ norm(SysMat'*DataVec);
    NoiseLevelCGLS = [x];
    %for i = 1:1:size(NoiseLevelCGLS,2)
    optionsHBLSQR.RegParam = 'discrep';%'WGCV';%
    %optionsHBLSQR.GCVweight = 'adapt';%0;%
    optionsHBLSQR.x_true = OrigVec;
    optionsHBLSQR.NoiseLevel = NoiseLevelCGLS;%0;
    optionsHBLSQR.eta = 1.01;
    optionsHBLSQR.RegMatrix = F;
    optionsHBLSQR.MaxIter = 200;
    %optionsHBLSQR.DecompOut  = 'on';
    optionsHBLSQR.Reorth = 'off';
    %dbstop in IRhybrid_lsqr
    tic
    [TikhImageHBLSQR, infoHBLSQR] = IRhybrid_lsqr(SysMat, DataVec, optionsHBLSQR);
    toc
    infoHBLSQR.TimeElaspsed = toc;
    TikhImageHBLSQRStore = TikhImageHBLSQR;
    %for j = 1:1:size(TikhImageHBLSQRStore,2)
    TikhImageHBLSQR = abs(reshape(TikhImageHBLSQRStore, ImgDim));
    ReconsImgHBLSQR = abs(TikhImageHBLSQR);
    ReconsImgHBLSQR = abs(reshape(ReconsImgHBLSQR, ImgDim));
    %figure(i), imagesc(ReconsImgHBLSQR)
    %pause(.01);
    %end
    %end
    %options  = IRhybrid_lsqr('defaults');
 
 %% Combining GCV and discrep &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
 %{
    GCVHBLSQRImage = TikhImageHBLSQRStore(:,size(TikhImageHBLSQRStore,2));
    GCVNoiseEst = norm(DataVec-(SysMat*GCVHBLSQRImage))/ norm(DataVec);
    optionsHBLSQR.RegParam = 'discrep';
    optionsHBLSQR.x_true = OrigVec;
    optionsHBLSQR.NoiseLevel =  GCVNoiseEst;
    optionsHBLSQR.eta = 1.01;
    optionsHBLSQR.RegMatrix = F;
    optionsHBLSQR.NE_Rtol = 10^(-17);
    optionsHBLSQR.MaxIter = 200;
    %dbstop in IRhybrid_lsqr
    [TikhImageHBLSQR, infoHBLSQR] = IRhybrid_lsqr(SysMat, DataVec, [1:200], optionsHBLSQR);
    TikhImageHBLSQRStore = TikhImageHBLSQR;
    for j = 1:1:size(TikhImageHBLSQR,2)
    TikhImageHBLSQR = abs(reshape(TikhImageHBLSQRStore(:,j), ImgDim));
    ReconsImg = abs(TikhImageHBLSQR);
    ReconsImg = abs(reshape(ReconsImg, ImgDim));
    imagesc(ReconsImg)
    pause(.01);
    end
    %}
 %{   
    % Xnrm Enrm Rnrm
    plot(infoHBLSQR.Xnrm,'LineWidth',3)
    %xlim([0 12])    
    hold on
    plot(infoHTV.Xnrm,'LineWidth',1.5)
    %ylim([0 1])
    hold off
    legend('HBLSQR','HTV')
 
    plot(b,'color',LColor{a},'linestyle',LStyle{a},'Marker',MStyle{a},'MarkerIndices',length(b),'MarkerSize',MSize{a},'LineWidth',1.5)
    set(gcf, 'Position',  [523,420,714,389])
    xlim([0 12])
    xlabel('No. of Iterations')
    ylabel('Relative Error Norm')
    title('IRhtv')
    hold on
    
Lap = F;
InvLap = pinv(F);
PIdnty = Lap\Lap;
I = eye(size(PIdnty));
TestOrg = OrigVec;
OprImg = TestOrg;
norm(TestOrg)
norm(OprImg)
for i=1:1:10
    OprImg = Lap*TestOrg;
    TestReconsImg = abs(reshape(OprImg, ImgDim));
    imagesc(TestReconsImg)
    pause(.01);
end

%}
    
 
end






