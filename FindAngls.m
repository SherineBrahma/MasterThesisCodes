close all;
clear all;

    %% ##################### ASSIGNING CONSTANTS ##############################
    
    global Gamma 

    Gamma = 2*pi*42.577e+06; % gyromagnetic ratio in rad/Ts

    %% ############## LOADING MEASURED MAGNETIC FIELD #########################

    [BZMsrd, XMsrdCoord, ZMsrdCoord] = LoadMsrdBZField("Bz_measured.mat");

    %% ########## QUADRATIC FITTING OF THE MEASURED MAGNETIC FIELD ############

    [BZFit, BZlsqCoeff] = QuadFitBZMsrd(BZMsrd, XMsrdCoord, ZMsrdCoord);

    %% ####### INTERPOLATING THE MEASURED AND FITTED MAGNETIC FIELDS ##########

    [BZFitIntrpol, BZMsrdIntrpol, XImgCoord, ZImgCoord, StpSizeXImg, StpSizeZImg] = IntrpolBZ(BZFit, BZMsrd, XMsrdCoord, ZMsrdCoord);
    %BZMsrdIntrpol = importdata("GradBMap.mat");  %Use this if you want to use some other B0 profile
    
    %% ######################## SET PARAMETERS ############################
     
    %Angle
    % Angle Domain
    AnglArryCell = {};
    AnglArryCell{1,1} = 0:5:355; % 0-360 Deg, 5 Deg Step
    % No of angle used
    NoOfAnglUsedArry = 36;
    
    %Time
    % Time step size
    DwTimeArry = 10*1e-06;       
    % No of time samples
    TDmnArry = 516;
    
    % No of Instances
    StatInstNo = 1;
    RecomAnglMat = [];
    
    % Load Gram matrix
    load('E:/G0To360Stp5');
    
for InsCnt = 1:1:StatInstNo
    RecomAnglMat = {};
    for AnglCnt = 1:1:size(NoOfAnglUsedArry,2)
        TDReconsImgCell = {};
        for TDCnt = 1:1:size(TDmnArry,2)
            DWReconsImgCell = {};
            for DWCnt = 1:1:size(DwTimeArry,2)
                %% ###################### LOOP Status #########################
                display(InsCnt)
                display(AnglCnt)
                display(TDCnt)
                display(DWCnt)
                %% ############# GETTING MATHEMATICAL MODEL OF SYSTEM #################
                
                % Angle sequences
                AnglArry = AnglArryCell{1,1};
                
                % Time step size
                StepSizeTime = DwTimeArry(DWCnt);        % Time step size is in seconds
                % No of time samples
                NoOfTimeInstances = TDmnArry(TDCnt);      % Total number of time instances considered
                
                disp('.... Getting mathematical model of the system....');

                % ******* SheppLoganMRINoise SheppLoganMRINoise50 SheppLoganMRINoNoise  data_total ********

                %Get parameters to generate system matrix and regularization matrix
                [AcqSMat, AcqPVectMapCell, AcqKVectMapCell, AcqKTpzmVectMapCell, AcqTimeArry, AcqStpSizeArea] = GetAcqSysMatCompParams( BZMsrd, BZMsrdIntrpol, XImgCoord, ZImgCoord, AnglArry, StepSizeTime, NoOfTimeInstances, StpSizeXImg, StpSizeZImg);

                %Get system matrix and regularization matrix
                [SysMat] = GetAcqSysMat(AcqPVectMapCell, AcqKVectMapCell, AcqTimeArry, AcqStpSizeArea, AnglArry);
            
                disp('.... Computing recommended angles....')
                
                NoOfAnglUsed = NoOfAnglUsedArry(AnglCnt);
                %G = GetNormGramMat(SysMat);
                RecomAnglArry = RecomAngls(G, NoOfAnglUsed, AnglArry, NoOfTimeInstances );
                RecomAnglMat{AnglCnt,1} = NoOfAnglUsed;
                RecomAnglMat{AnglCnt,2} = RecomAnglArry;                
            end
        end
    end
end
     
    
