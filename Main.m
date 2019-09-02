    close all;
    clear all;

    %% ##################### ASSIGNING CONSTANTS ##############################
    
    global Gamma 

    Gamma = 2*pi*42.577e+06; % gyromagnetic ratio in rad/Ts

    %% ############## LOADING MEASURED MAGNETIC FIELD #########################

    [BZMsrd, XMsrdCoord, ZMsrdCoord] = LoadMsrdBZField("Bz_measured.mat");

    %% ########## QUADRATIC FITTING OF THE MEASURED MAGNETIC FIELD ############

    [BZFit, BZlsqCoeff] = QuadFitBZMsrd(BZMsrd, XMsrdCoord, ZMsrdCoord);

    %% ############### MEASURED VS FITTED MAGNETIC FLUX PLOT ##################
    %{
    figure(1);
    subplot(1,2,1);imagesc( XMsrdCoord, ZMsrdCoord, BZMsrd);
    axis xy;
    colorbar;
    title('Measured z-component magnetic flux density');
    subplot(1,2,2);
    imagesc( XMsrdCoord, ZMsrdCoord, BZFit);
    axis xy;
    colorbar;
    title('Fitted z-component magnetic flux density');
    %}
    %% ####### INTERPOLATING THE MEASURED AND FITTED MAGNETIC FIELDS ##########

    [BZFitIntrpol, BZMsrdIntrpol, XImgCoord, ZImgCoord, StpSizeXImg, StpSizeZImg] = IntrpolBZ(BZFit, BZMsrd, XMsrdCoord, ZMsrdCoord);
    
    %% ##### MEASURED VS FITTED INTERPOLATED MAGNETIC FLUX PLOT ###############
    %{
    figure(2);
    subplot(1,2,1);
    imagesc( XImgCoord, ZImgCoord, BZMsrdIntrpol);
    axis xy;
    colorbar;
    title('Interpolated z-component magnetic flux density')
    subplot(1,2,2);
    imagesc(XImgCoord, ZImgCoord, BZFitIntrpol);
    axis xy;
    colorbar;
    title('Interpolated fitted z-component magnetic flux density')
    %}
    
    %% ############## FIELD VARIATIONS USING QUADRATIC FIT ################
    %{
    GradBZFit = GetGradBZFit(BZlsqCoeff, XImgCoord, ZImgCoord, StpSizeXImg, StpSizeZImg);
    figure(3);
    imagesc(GradBZFit);
    colorbar;
    %}
    
    %% ######################## SET PARAMETERS ############################
     
    %Angle
    AnglArryCell = {};
    AnglArryCell{1,1} = 0:10:350; % 0-360 Deg, 10 Deg Step
    AnglArryCell{2,1} = 0:5:355; % 0-360 Deg, 5 Deg Step
    AnglArryCell{3,1} = 0:10:170; % 0-180 Deg, 10 Deg Step
    AnglArryCell{4,1} = 0:5:175; % 0-180 Deg, 5 Deg Step
    AnglArryCell{5,1} = 0:2.5:177.5; % 0-180 Deg, 2.5 Deg Step
    AnglArryCell{6,1} = 0:10:80; % 0-90 Deg, 10 Deg Step
    AnglArryCell{7,1} = 0:5:85; % 0-90 Deg, 5 Deg Step
    AnglArryCell{8,1} = 0:2.5:87.5; % 0-90 Deg, 2.5 Deg Step
    
    %Time
    StepSizeTime = 10e-06;        % Time step size is in seconds
    NoOfTimeInstances = 516;      % Total number of time instances considered
   
    ExpNo = size(AnglArryCell, 1);    
    
StatInstNo = 20;
for k = 1:1:StatInstNo
    display(k)
    ReconsImgCell = {};
for i = 1:1:ExpNo
    display(i)
    %% ###################### LOOP INITIALIZATION #########################
    
    clear DataStore;
    
    %% ############# GETTING MATHEMATICAL MODEL OF SYSTEM #################
    
    AnglArry = AnglArryCell{i,1};
    
    disp('.... Getting mathematical model of the system....');
        
    % ******* SheppLoganMRINoise SheppLoganMRINoise50 SheppLoganMRINoNoise  data_total ********
    
    %Get parameters to generate system matrix and regularization matrix
    [AcqSMat, AcqPVectMapCell, AcqKVectMapCell, AcqKTpzmVectMapCell, AcqTimeArry, AcqStpSizeArea] = GetAcqSysMatCompParams( BZMsrd, BZMsrdIntrpol, XImgCoord, ZImgCoord, AnglArry, StepSizeTime, NoOfTimeInstances, StpSizeXImg, StpSizeZImg);
    
    %Get system matrix and regularization matrix
    [AcqSysMat] = GetAcqSysMat(AcqPVectMapCell, AcqKVectMapCell, AcqTimeArry, AcqStpSizeArea, AnglArry);
  


    %% #################### Simulate IMAGE #################################
    
    disp('.... Preparing sample....');
    m = size(BZFitIntrpol,1);
    n = m;
    Img = phantom('Shepp-Logan',m);                                % ## Inbuilt function to create shepp-logan phantom
    Img(Img>=0.1)=0.04;                                            % ## Some thresholding done
    Img = imresize(Img,[m n]);                                     % ## Does some resizing
    ImgStkd  = Img(:);                                             % ## Stacks the image as a column vector
    disp('.... Imaging Sample....');
    SimData = AcqSysMat*ImgStkd;                                   % ## Measured data
    SNR = 10;
    SimData = awgn(SimData, SNR,'measured','dB');                  % Adding Noise in dB scale 'dB'
    
    %{
    %% Test to see lsqr results
    ImgStkd = lsqr(SysMat, SimData);
    imagesc(abs(reshape(ImgStkd,65,65)))
    %}    
    
    %{
    %% Test to see whether SNR function is working
    ExpNo = 10000;
    for i = 1:1:ExpNo
    A = [1 1 1 1];
    N = [1 1 1 1];
    SNRCal = 10 * log10(norm(A)/norm(N));
    SNR = 0;
    TestData = awgn(A, SNR,'measured','dB');
    EstNoise = TestData - A;
    EstNoiseNorm(i) = norm(EstNoise);
    end
    hist(EstNoiseNorm, ExpNo/100)
    %}
    
    %{
    %% Test to see whether the signals are decaying
    AngleNo = 7;
    LBndSamp = 1 + NoOfTimeInstances*(AngleNo-1);
    UBndSamp = NoOfTimeInstances*(AngleNo-1) + NoOfTimeInstances;
    SimDataTest = SimData( LBndSamp: UBndSamp);
    figure, plot(abs(SimDataTest))
    %figure, plot(real(SimDataTest))
    %figure, plot(imag(SimDataTest))
    %}
    
    SimDataUnstkd = reshape(SimData,NoOfTimeInstances,[]);
    RealSimDataUnstkd = real(SimDataUnstkd);
    ImgSimDataUnstkd = imag(SimDataUnstkd);
    DataStore(:,1,:) = RealSimDataUnstkd;
    DataStore(:,2,:) = ImgSimDataUnstkd;

    %% #################### GET DATA CELL  ###########################
    
    %Get parameters to generate system matrix and regularization matrix
    [DataCell] = GetDataMatCompParams(DataStore, AnglArry);
    
    %% #################### GET COMPUTATION MATRICES ###########################
    
    [ SMat, PVectMapCell, KVectMapCell, KTpzmVectMapCell, TimeArry, StpSizeArea] = GetCompMatParams( BZMsrd, BZMsrdIntrpol, XImgCoord, ZImgCoord, AnglArry, StepSizeTime, NoOfTimeInstances, StpSizeXImg, StpSizeZImg);
    disp('.... Getting Computation Matrices....');
    %Get matrices needed for computations
    [SysMat, DataVec, F] = GetCompMat(PVectMapCell, KVectMapCell, TimeArry, StpSizeArea, AnglArry, DataCell);
    
    %% #################### COMPUTE IMAGE #################################
    
    disp('.... Computing image....');
    %{
    % Reconstruct using pseudo inverse of the system matrix
    ReconsImg = PInvRecon(PVectMapCell, KVectMapCell, DataCell, TimeArry, StpSizeArea);
    imagesc(ReconsImg)
    title('Using Pseudo Inverse')
    
    % Reconstruct using adjoint operator matrix  
    ReconsImg = AdjOpRecon(PVectMapCell, KTpzmVectMapCell, DataCell, SMat, TimeArry, AnglArry, StepSizeTime, StepSizeAngl);
    figure, imagesc(ReconsImg)
    title('Using Adjoint Operator')
    
    % Reconstruct using least square minimization
    ReconsImg = LsqrRecon(PVectMapCell, KVectMapCell, DataCell, TimeArry, StpSizeArea);
    figure, imagesc(ReconsImg)
    title('Using Least Square Minimization')
        
    % Reconstruct using Tikhonov Regularization
    Lambda = .1;
    ReconsImg = TikhRecon(Lambda, PVectMapCell, KVectMapCell, DataCell, TimeArry, StpSizeArea);
    figure, imagesc(ReconsImg)
    title('Using Tikhonov Regularization')
    
    % Reconstruct using Total Variation Regularization
    Lambda = .1;
    ReconsImg = TotalVarRecon(Lambda, PVectMapCell, KVectMapCell, DataCell, TimeArry, StpSizeArea);
    figure, imagesc(ReconsImg)
    title('Using Tikhonov Regularization')
    
    % Reconstructing using MLE is equivalent to least square minimization
    
    % Reconstruct using Huber-Markov Random Field Model Regularization
    Lambda = .1;
    ReconsImg = HuberMarkovRecon(Lambda, PVectMapCell, KVectMapCell, DataCell, TimeArry, StpSizeArea);
    figure, imagesc(ReconsImg)
    title('Using Tikhonov Regularization')
    %}
    
    % Reconstruction using IR Tools
    Lambda = .1;
    [ReconsImgHTV, infoHTV, ReconsImgHBLSQR, infoHBLSQR] = IRTikhRecon(SysMat, F, DataVec, Img);
    %figure, imagesc(ReconsImg)
    %title('Using Tikhonov Regularization')
    
    ReconsImgCell{i,1} = ReconsImgHBLSQR;
    ReconsImgCell{i,2} = infoHBLSQR;
    ReconsImgCell{i,3} = ReconsImgHTV;
    ReconsImgCell{i,4} = infoHTV;
end
   StatReconsImgCell{k,1} = ReconsImgCell;
end
save('AngleInvestigation', 'StatReconsImgCell')
%{

InsNo = 6;
LStyle = {'-','-','-','-','-','-','-','--'};
for InsNo=1:1:20
AngoNo = 2;
for g = 1:1:ExpNo
    plot(StatReconsImgCell{InsNo,1}{g,2*AngoNo}.Enrm,'linestyle',LStyle{g},'LineWidth',2)
    xlim([0 10])
    hold on
end
hold off
%legend('0-180 Deg, 5 Deg Step','0-180 Deg, 2.5 Deg Step','0-90 Deg, 10 Deg Step','0-90 Deg, 5 Deg Step','0-90 Deg, 2.5 Deg Step')
legend('0-360 Deg, 10 Deg Step','0-360 Deg, 5 Deg Step','0-180 Deg, 10 Deg Step','0-180 Deg, 5 Deg Step','0-180 Deg, 2.5 Deg Step','0-90 Deg, 10 Deg Step','0-90 Deg, 5 Deg Step','0-90 Deg, 2.5 Deg Step')
pause;
end

%}

%{
SumReconsImgMat = zeros(size(ReconsImgCell{1,1}));
for i = 1:1:ExpNo
    
    IntSumReconsImgCell = ReconsImgCell{i,1};
    SumReconsImgMat = SumReconsImgMat + IntSumReconsImgCell;
    
end
SumReconsImgMat = SumReconsImgMat/ExpNo;
figure, imagesc(SumReconsImgMat)
%}     
    
