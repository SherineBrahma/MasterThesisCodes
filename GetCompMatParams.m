function [ SMat, PVectMapCell, KVectMapCell, KTpzmVectMapCell, TimeArry, StpSizeArea] = GetCompMatParams( BZMsrd, BZMsrdIntrpol, XImgCoord, ZImgCoord, AnglArry, StepSizeTime, NoOfTimeInstances, StpSizeXImg, StpSizeZImg)

    global Gamma 
    
    %% Getting Area Step Size
    StpSizeArea = StpSizeXImg * StpSizeZImg;
    
    %% Making TimeArry and AnglArry
    % Getting the time instances used for computations
    TimeArry = (0 : StepSizeTime: StepSizeTime*(NoOfTimeInstances-1))';    
    % Getting the angles used for computations
    AnglArry = AnglArry.';
    NoOfAnglInstances = size(AnglArry,1);
    AngleCell = num2cell( AnglArry, 2);
    
    %% Delta omega Map
    % Making coordinates Grids
    [ XImgMesh, ZImgMesh] = meshgrid( XImgCoord, ZImgCoord);
    YImgMesh = zeros(size(ZImgMesh)); % To be used later during creating magnetization map
    RotMeshCell = {};
    RotMeshCell{1,1} = XImgMesh;
    RotMeshCell{1,2} = ZImgMesh;
    RotMeshCell = repmat(RotMeshCell, NoOfAnglInstances,1);
    % Rotating the Mesh Grids of x and z imaging coordinates by corresponding angles
    RotMeshCell( :, 3) = cellfun(@(x,z, theta)  (x * cosd(theta)) - (z * sind(theta)), RotMeshCell(:,1) , RotMeshCell(:,2), AngleCell, 'UniformOutput',false);
    RotMeshCell( :, 4) = cellfun(@(x,z, theta)  (x * sind(theta)) + (z * cosd(theta)), RotMeshCell(:,1) , RotMeshCell(:,2), AngleCell, 'UniformOutput',false);
    % Getting the magnetic fields in the obtained rotated meshes
    CompOmegaMapCell( :, 1) = cellfun(@(x,z)  interp2( XImgMesh, ZImgMesh, BZMsrdIntrpol, x, z, 'spline'), RotMeshCell(:,3) , RotMeshCell(:,4), 'UniformOutput',false);
    % Getting delta omega in the obtained rotated meshes
    BZMsrdCenter = BZMsrdIntrpol(floor((size(BZMsrdIntrpol,1))/2) + 1, floor((size(BZMsrdIntrpol,2))/2) + 1); % Getting the field strength of the central coordinates of the measured field
    DelOmegaMapCell = cellfun(@(B)  Gamma * (B - BZMsrdCenter) , CompOmegaMapCell(:,1), 'UniformOutput',false);
    %save('C:\Users\Sherine\Desktop\DelftStudy\Thesis\Codes\DataSimulationCode\DelOmegaMapCell','DelOmegaMapCell');
    %rankcell = cellfun(@(SysMatAngl)  rank(SysMatAngl) , DelOmegaMapCell, 'UniformOutput',false);
   %% T2 Map
    % Getting the T2 Map
    T2 = 90e+09;
    T2MapCell = cellfun(@(B)  T2 * ones(size(B)) , DelOmegaMapCell, 'UniformOutput',false);
    
    %Experiment
    %A = imread("BlurCircle.png");
    %A = rgb2gray(A);
    %A = double(A);
    %A = imresize(A,[65 65]);
    %A = A./(max(max(A))-min(min(A)));
    %imagesc(A)
    %mean(mean(A))
    
   %% Magnetization Map
    
    ImgDim = 65;
    MapMsk = ones(size(BZMsrdIntrpol));
    B1Strength = 2 * 10^(-6); %In Tesla
    Omega1 = Gamma*B1Strength;
    Mz0Equib = 1; 
    
    % Magnetization after flip by alpha pulse
    Alpha = 45;
    Mx0Map = 0 * MapMsk;
    My0Map = Mz0Equib * sind(Alpha) * MapMsk;
    Mz0Map = Mz0Equib * cosd(Alpha) * MapMsk;
    
    % Adding noise to the flipped angles
    
    AngleNoise = [0 0];  % No angle noise is used during reconstruction because we do not know it 
    
    Mx0MapStkd = reshape(Mx0Map,1,[]);
    My0MapStkd = reshape(My0Map,1,[]);
    Mz0MapStkd = reshape(Mz0Map,1,[]);
    M0MapStkd = [Mx0MapStkd; My0MapStkd; Mz0MapStkd];
    
    CompM0NoiseCell = repmat({M0MapStkd},[NoOfAnglInstances 1]);
    CompM0NoiseCell( :, 2) = cellfun(@(M0)  AddNoiseM0Map( M0, AngleNoise), CompM0NoiseCell( :, 1), 'UniformOutput',false);
    M0NoiseCell( :, 1) = cellfun(@(M0)  reshape(M0(1,:),ImgDim,ImgDim), CompM0NoiseCell( :, 2), 'UniformOutput',false);
    M0NoiseCell( :, 2) = cellfun(@(M0)  reshape(M0(2,:),ImgDim,ImgDim), CompM0NoiseCell( :, 2), 'UniformOutput',false);
    M0NoiseCell( :, 3) = cellfun(@(M0)  reshape(M0(3,:),ImgDim,ImgDim), CompM0NoiseCell( :, 2), 'UniformOutput',false);
    
    % Initial phase due to rotating coordinates to room coordinates is ignored
    
    % Magnetization For Rotated Setting
    RotdMagMap( :, 1) = cellfun(@(Mx, My, Mz, theta)  (Mx * cosd(theta)) - (Mz * sind(theta)), M0NoiseCell( :, 1), M0NoiseCell( :, 2), M0NoiseCell( :, 3), AngleCell, 'UniformOutput',false);
    RotdMagMap( :, 2) = cellfun(@(Mx, My, Mz, theta)  My, M0NoiseCell( :, 1), M0NoiseCell( :, 2), M0NoiseCell( :, 3), AngleCell, 'UniformOutput',false);
    RotdMagMap( :, 3) = cellfun(@(Mx, My, Mz, theta)  (Mx * sind(theta)) + (Mz * cosd(theta)), M0NoiseCell( :, 1), M0NoiseCell( :, 2), M0NoiseCell( :, 3), AngleCell, 'UniformOutput',false);
    
    %{
    % Magnetization Visualization
    AngNoMag = 1;
    ViewMxMap = RotdMagMap{AngNoMag,1};%M0NoiseCell{ 1, 1};%Mx0Map;%
    ViewMyMap = RotdMagMap{AngNoMag,2};%M0NoiseCell{ 1, 2};%My0Map;%
    ViewMzMap = RotdMagMap{AngNoMag,3};%M0NoiseCell{ 1, 3};%Mz0Map;%
    Step = 4;
    figure, quiver3(XImgMesh(1:Step:ImgDim,1:Step:ImgDim), YImgMesh(1:Step:ImgDim,1:Step:ImgDim), -ZImgMesh(1:Step:ImgDim,1:Step:ImgDim), ViewMxMap(1:Step:ImgDim,1:Step:ImgDim), ViewMyMap(1:Step:ImgDim,1:Step:ImgDim), ViewMzMap(1:Step:ImgDim,1:Step:ImgDim),0.5)
    xlabel('x-axis')
    ylabel('y-axis')
    zlabel('z-axis')
    axis([-.055 .055 -.055 .055 -0.055 .055])
    %}
    
    
    % Getting the Magnetization Map( M = Mxcos(theta) - Mzsin(theta) + iMy)
    CompMagtznMapCell = RotdMagMap;
    MagtznMapCell= cellfun(@( x, y, z, theta)  (x * cosd(theta)) - (z * sind(theta)) + (1i * y), CompMagtznMapCell(:, 1) , CompMagtznMapCell(:, 2), CompMagtznMapCell(:, 3), AngleCell, 'UniformOutput',false);
   
   %% Receive Field Map
    % Receive Field x-component Map
    WxMap = 0 * ones(size(BZMsrdIntrpol));
    CompRcvFldMapCell( :, 1) = cellfun(@(x,z)  interp2( XImgMesh, ZImgMesh, WxMap, x, z, 'spline'), RotMeshCell(:,3) , RotMeshCell(:,4), 'UniformOutput',false);
    % Receive Field y-component Map
    WyMap = 1 * ones(size(BZMsrdIntrpol));
    CompRcvFldMapCell( :, 2) = cellfun(@(x,z)  interp2( XImgMesh, ZImgMesh, WyMap, x, z, 'spline'), RotMeshCell(:,3) , RotMeshCell(:,4), 'UniformOutput',false);
    % Receive Field z-component Map
    WzMap = 0 * ones(size(BZMsrdIntrpol));
    CompRcvFldMapCell( :, 3) = cellfun(@(x,z)  interp2( XImgMesh, ZImgMesh, WzMap, x, z, 'spline'), RotMeshCell(:,3) , RotMeshCell(:,4), 'UniformOutput',false);
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
    
end