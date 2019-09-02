function [ DataCell, SMat, PVectMapCell, KVectMapCell, KTpzmVectMapCell, TimeArry, AnglArry, StpSizeArea] = GetImageCompParams(MRIDataStoreFile, BZMsrd, BZMsrdIntrpol, XImgCoord, ZImgCoord, StepSizeAngl, NoOfAnglInstances, StepSizeTime, NoOfTimeInstances, StpSizeXImg, StpSizeZImg)

global Gamma 

   %% Reading Data

    MRIDataStore = load(MRIDataStoreFile);             % Load Dataset

    %Chopping the data sequences to include only those after the peak

    % Format of DataVect:
    % Row -> Angles Axis
    % Column -> Time Axis
    RawDataVect = (squeeze(MRIDataStore.data_store(:,1,:)) + (1i * squeeze(MRIDataStore.data_store(:,2,:)))).'; 
    %plot(abs(RawDataVect(5,:)))

    %% Getting Area Step Size
    
    StpSizeArea = StpSizeXImg * StpSizeZImg;
    
    %% Making TimeArry and AnglArry
    % Finding the peak of w.r.t. time axis
    %[~ , PeakIndxArry] = max(abs( RawDataVect), [], 2);
    %PeakIndx = max(PeakIndxArry);
    %PeakIndx = min(PeakIndx, (size(RawDataVect, 2) - NoOfTimeInstances)); % This is to handle when PeakIndx is close the end of the array index
    %TimeIndxArry = PeakIndx : (PeakIndx + NoOfTimeInstances - 1);
    %TimeArry = TimeIndxArry * StepSizeTime;
    %TimeArry = TimeArry - TimeArry( 1, 1);

    TimeArry = MRIDataStore.time_mu_s.';
    
    % Getting the angles used for acquisition
    
    AnglArry = 0: StepSizeAngl: (StepSizeAngl * (NoOfAnglInstances - 1));
    AnglArry = AnglArry.';
    AngleCell = num2cell( AnglArry, 2);
    
    %AnglArry = AnglArry.*(2*pi/360);      % Some Processing of angles which I am not sure of
    %AnglArry = AnglArry + pi/2;
    %AnglArry = -AnglArry;
    %AnglArry = AnglArry.*(360/(2*pi));
    % Truncating the RawDataVect
    %DataVect = RawDataVect( 1: size( AnglArry, 1), TimeIndxArry );
    
   %% DataCell
    % Getting DataCell
    DataVect = RawDataVect( 1: size( AnglArry, 1), :);
    DataCell = num2cell(DataVect, 2);
    DataCell = cellfun(@(x)  x.', DataCell, 'UniformOutput',false);       % Putting the same time varying data to be processed for different imaging angle rotating matrix
   
   %% Delta omega Map
    % Making coordinates Grids
    [ XImgMesh, ZImgMesh] = meshgrid( XImgCoord, ZImgCoord);    
    CompOmegaMapCell = {};
    CompOmegaMapCell{1,1} = XImgMesh;
    CompOmegaMapCell{1,2} = ZImgMesh;
    CompOmegaMapCell = repmat(CompOmegaMapCell, NoOfAnglInstances,1);
    % Rotating the Mesh Grids of x and z imaging coordinates by corresponding angles
    CompOmegaMapCell( :, 3) = cellfun(@(x,z, theta)  (x * cosd(theta)) - (z * sind(theta)), CompOmegaMapCell(:,1) , CompOmegaMapCell(:,2), AngleCell, 'UniformOutput',false);
    CompOmegaMapCell( :, 4) = cellfun(@(x,z, theta)  (x * sind(theta)) + (z * cosd(theta)), CompOmegaMapCell(:,1) , CompOmegaMapCell(:,2), AngleCell, 'UniformOutput',false);
    % Getting the magnetic fields in the obtained rotated meshes
    CompOmegaMapCell( :, 5) = cellfun(@(x,z)  interp2( XImgMesh, ZImgMesh, BZMsrdIntrpol, x, z, 'spline'), CompOmegaMapCell(:,3) , CompOmegaMapCell(:,4), 'UniformOutput',false);
    % Getting delta omega in the obtained rotated meshes
    BZMsrdCenter = BZMsrd(floor((size(BZMsrd,1))/2) + 1, floor((size(BZMsrd,2))/2) + 1);     % Getting the field strength of the central coordinates of the measured field
    DelOmegaMapCell = cellfun(@(B)  Gamma * (B - BZMsrdCenter) , CompOmegaMapCell(:,5), 'UniformOutput',false);
    %save('C:\Users\Sherine\Desktop\DelftStudy\Thesis\Codes\DataSimulationCode\DelOmegaMapCell','DelOmegaMapCell');
  
   %% Flip Angle Map
    Alpha = 90;
    Tp = (Alpha * pi)/(180 * Gamma * BZMsrdCenter);
    FlipAngMapCell( :, 1) = cellfun(@(B0)  (180 * Gamma * B0 * Tp)/pi, CompB0MapCell( :, 1), 'UniformOutput',false);
    
   %% T2 Map
    % Getting the T2 Map
    T2 = 90e+09;
    T2MapCell = cellfun(@(B)  T2 * ones(size(B)) , DelOmegaMapCell, 'UniformOutput',false);
    
   %% Magnetization Map
   
    ImgDim = 65;
    M0MapMsk = ones(size(BZMsrdIntrpol));
    Mx0Map = 0 * M0MapMsk;
    My0Map = 0 * M0MapMsk;
    Mz0Map = 1 * M0MapMsk;
    
    % Equilibrium Magnetization For Rotated Setting
    RotdMagMap( :, 1) = cellfun(@(theta)  (Mx0Map * cosd(theta)) - (Mz0Map * sind(theta)), AngleCell, 'UniformOutput',false);
    RotdMagMap( :, 2) = cellfun(@(theta)  My0Map,  AngleCell, 'UniformOutput',false);
    RotdMagMap( :, 3) = cellfun(@(theta)  (Mx0Map * sind(theta)) + (Mz0Map * cosd(theta)), AngleCell, 'UniformOutput',false);
    
    % Magnetization Visualization
    %{
    AngNoMag = 1;
    ViewMxMap = RotdMagMap{AngNoMag,1};
    ViewMyMap = RotdMagMap{AngNoMag,2};
    ViewMzMap = RotdMagMap{AngNoMag,3};
    Step = 4;
    quiver3(XImgMesh(1:Step:ImgDim,1:Step:ImgDim), YImgMesh(1:Step:ImgDim,1:Step:ImgDim), -ZImgMesh(1:Step:ImgDim,1:Step:ImgDim), ViewMxMap(1:Step:ImgDim,1:Step:ImgDim), ViewMyMap(1:Step:ImgDim,1:Step:ImgDim), ViewMzMap(1:Step:ImgDim,1:Step:ImgDim),0.5)
    xlabel('x-axis')
    ylabel('y-axis')
    zlabel('z-axis')
    axis([-.055 .055 -.055 .055 -0.055 .055])
    %}
    % Maps after flipping about the x' axis
    
    %(Left to do)
    
    % Getting the Magnetization Map( M = Mxcos(theta) - Mzsin(theta) + iMy)
    CompMagtznMapCell = RotdMagMap; % Quick fix. Mend it later.
    MagtznMapCell= cellfun(@( x, y, z, theta)  (x * cosd(theta)) - (z * sind(theta)) + (1i * y), CompMagtznMapCell(:, 1) , CompMagtznMapCell(:, 2), CompMagtznMapCell(:, 3), AngleCell, 'UniformOutput',false);
   
   %% Receive Field Map
    % Receive Field x-component Map
    Wx = 0;
    CompRcvFldMapCell( :, 1) = cellfun(@(B)  Wx * ones(size(B)) , DelOmegaMapCell, 'UniformOutput',false);
    % Receive Field y-component Map
    Wy = 1i;
    CompRcvFldMapCell( :, 2) = cellfun(@(B)  Wy * ones(size(B)) , DelOmegaMapCell, 'UniformOutput',false);
    % Receive Field z-component Map
    Wz = 0;
    CompRcvFldMapCell( :, 3) = cellfun(@(B)  Wz * ones(size(B)) , DelOmegaMapCell, 'UniformOutput',false);
    % Getting the Receive Field (  W = Wxcos(theta) - Wzsin(theta) - iWy  )
    RcvFldMapCell = cellfun(@( x, y, z, theta)  (x * cosd(theta)) - (z * sind(theta)) - (1i * y), CompRcvFldMapCell(:, 1) , CompRcvFldMapCell(:, 2), CompRcvFldMapCell(:, 3), AngleCell, 'UniformOutput',false);
    
   %% P Vector Map
    % Getting the P Vector Map
    PVectMapCell = cellfun(@( T2, DelOmega) (1./T2) + (1i*DelOmega), T2MapCell, DelOmegaMapCell, 'UniformOutput',false);
    
   %% K Vector Map
    % Getting the K Vector Map(The constants are ignored)
    KVectMapCell = cellfun(@( M, W)  -1i * (M .* W), MagtznMapCell , RcvFldMapCell, 'UniformOutput',false);
    
    % Getting the K' vector (Used for employing Trapezium rule)
    KTpzmVectMapCell = KVectMapCell;
    KTpzmVectMapCell{ 1, 1} = KTpzmVectMapCell{ 1, 1}./2;
    KTpzmVectMapCell{ size(KTpzmVectMapCell,1), 1} = KTpzmVectMapCell{ size(KTpzmVectMapCell,1), 1}./2;
   
   %% SMat
    %Getting SMat
    SMat = ones(size(TimeArry'));
    SMat(1,:) = 1/2;
    SMat(size(TimeArry',1),:) = 1/2;
    SMat = diag(SMat);    
    
    %{
    IntrCompParamsCell = cellfun(@(data) SMat * exp(-( * TimeArry')) * data, DataCell, 'UniformOutput',false);
    CompSysMatCell2 = {};
    CompSysMatCell2 = CompParamsMapCell( :, 15);                    % p
    CompSysMatCell2( :, 2) = CompParamsMapCell( :, 19);             % K'
    CompSysMatCell2( :, 3) = cellfun(@(x)  num2cell(x), CompSysMatCell2( :, 1), 'UniformOutput',false);
    CompSysMatCell2( :, 4) = cellfun(@(x)  num2cell(x), CompSysMatCell2( :, 2), 'UniformOutput',false);
    CompSysMatCell2( :, 5) = cellfun(@(i, j) cellfun(@(k, l) (StepSizeAngl * StepSizeTime * SMat * l * exp(-(k * TimeArry'))), i, j, 'UniformOutput',false) , CompSysMatCell2( :, 3), CompSysMatCell2( :, 4), 'UniformOutput',false);

    % Getting -p * t
    % Format of PVect and KVect:
    % Row -> Angles Axis
    % Column -> Space Axis
    SysCell = CompSysMatCell2( :, 5);
    %}

end