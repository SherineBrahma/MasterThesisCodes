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
    %BZMsrdIntrpol = BZFitIntrpol;%importdata("GradBMap.mat");  %Use this if you want to use some other B0 profile
    
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
    %load('SuggestedAngles')
    AnglArryOptns = 0:5:355; % 0-360 Deg, 5 Deg Step
    AnglArryCell = {};
    AnglArryCell{1,1} = 0:10:350;
    AnglArryCell{2,1} = 0:5:355;
    %.6 Repeated
    %AnglArryCell{3,1} = [55,85,95,125,210,235,265,275,355,55,85,95,125,210,235,265,275,355,55,85,95,125,210,235,265,275,355,55,85,95,125,210,235,265,275,355];
    %.6 Extended
    %AnglArryCell{4,1} = [5,15,30,35,50,55,60,75,85,95,105,120,125,140,145,150,165,175,185,195,210,215,230,235,240,255,265,275,285,300,305,320,325,330,345,355];
    %AnglArryCell{1,1} = [5,125,145,215,235];
    %AnglArryCell{2,1} = [5,95,125,145,215,235,265,285,310,355];
    %AnglArryCell{3,1} = [5,15,30,40,55,75,80,100,105,120,125,130,150,160,175,195,230,250,265,275,285,320,325,350,355];
    %AnglArryCell{4,1} = [5,20,30,40,50,55,75,85,95,120,130,135,140,150,165,175,185,190,195,205,215,220,240,245,250,255,260,280,285,300,305,310,320,325,340,355];
    %AnglArryCell{5,1} = [5,10,15,35,40,50,55,65,75,80,85,100,105,115,120,130,140,145,150,155,160,165,170,175,190,195,200,210,215,230,235,240,250,260,265,275,280,290,300,305,320,330,345,350,355];
    %AnglArryCell{6,1} = [5,10,15,25,30,35,45,50,55,60,70,80,85,95,100,105,110,115,120,125,130,140,145,150,155,165,170,175,185,190,195,200,210,215,220,230,240,245,255,260,265,275,280,285,295,300,305,320,325,330,335,340,345,350,355];
    %AnglArryCell{7,1} = [5,10,15,20,25,30,35,40,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,140,145,150,155,160,165,170,175,185,190,195,200,205,210,215,220,230,235,245,250,255,260,265,275,280,285,290,295,300,310,315,320,325,330,335,340,345,350,355];
    %AnglArryCell{8,1} = [0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150,155,160,165,170,175,180,185,190,195,200,205,210,215,220,225,230,235,240,245,250,255,260,265,270,275,280,285,290,295,300,305,310,315,320,325,330,335,340,345,350,355];
    %RndPerm = randperm(72);
    %AnglArryCell{4,1} = sort(AnglArryOptns(RndPerm(1:36)));
        
    %
    %Time
    % Time step size
    DwTimeArry = [10 1]*1e-06;        
    % No of time samples
    TDmnArry = [512 256];
    
    % No of Instances
    StatInstNo = 20;
    
for InsCnt = 1:1:StatInstNo
    AReconsImgCell = {};
    for AnglCnt = 1:1:size(AnglArryCell,1)
        TDReconsImgCell = {};
        for TDCnt = 1:1:size(TDmnArry,2)
            DWReconsImgCell = {};
            for DWCnt = 1:1:size(DwTimeArry,2)
                display(InsCnt)
                display(AnglCnt)
                display(TDCnt)
                display(DWCnt)
                %% ###################### LOOP INITIALIZATION #########################

                clear DataStore;

                %% ############# GETTING MATHEMATICAL MODEL OF SYSTEM #################

                AnglArry = AnglArryCell{AnglCnt,1};%Time
                % Time step size
                StepSizeTime = DwTimeArry(DWCnt);        % Time step size is in seconds
                % No of time samples
                NoOfTimeInstances = TDmnArry(TDCnt);      % Total number of time instances considered
                
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
                [ReconsImgHTV, infoHTV, ReconsImgHBLSQR, infoHBLSQR] = IRTikhRecon(SysMat, F, DataVec, Img);
                %figure, imagesc(ReconsImg)
                %title('Using Tikhonov Regularization')

                DWReconsImgCell{DWCnt,1} = StepSizeTime;
                DWReconsImgCell{DWCnt,2} = ReconsImgHBLSQR;
                DWReconsImgCell{DWCnt,3} = infoHBLSQR;
                DWReconsImgCell{DWCnt,4} = ReconsImgHTV;
                DWReconsImgCell{DWCnt,5} = infoHTV;
                
            end
            TDReconsImgCell{TDCnt,1} = NoOfTimeInstances;
            TDReconsImgCell{TDCnt,2} = DWReconsImgCell;
        end
        AReconsImgCell{AnglCnt,1} = AnglArry;
        AReconsImgCell{AnglCnt,2} = TDReconsImgCell;
    end
    StatReconsImgCell{InsCnt,1} = AReconsImgCell;
end
save('PSelectDiscrepInvest10dB', 'StatReconsImgCell')                
%save('AngleInvestigation', 'StatReconsImgCell')

%{
SumReconsImgMat = zeros(size(ReconsImgCell{1,1}));
for i = 1:1:ExpNo
    
    IntSumReconsImgCell = ReconsImgCell{i,1};
    SumReconsImgMat = SumReconsImgMat + IntSumReconsImgCell;
    
end
SumReconsImgMat = SumReconsImgMat/ExpNo;
figure, imagesc(SumReconsImgMat)
%}     
