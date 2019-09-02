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
    BZMsrdCenter = BZMsrd(floor((size(BZMsrd,1))/2) + 1, floor((size(BZMsrd,2))/2) + 1);     % Getting the field strength of the central coordinates of the measured field
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
    % Initial Magnetization x-component Map
    Mx0Map = 1 * ones(size(BZMsrdIntrpol));
    CompMagtznMapCell( :, 1) = cellfun(@(x,z)  interp2( XImgMesh, ZImgMesh, Mx0Map, x, z, 'spline'), RotMeshCell(:,3) , RotMeshCell(:,4), 'UniformOutput',false);
    % Initial Magnetization y-component Map
    My0Map = 0 * ones(size(BZMsrdIntrpol));
    CompMagtznMapCell( :, 2) = cellfun(@(x,z)  interp2( XImgMesh, ZImgMesh, My0Map, x, z, 'spline'), RotMeshCell(:,3) , RotMeshCell(:,4), 'UniformOutput',false);
    % Initial Magnetization z-component Map
    Mz0Map = 1 * ones(size(BZMsrdIntrpol));
    CompMagtznMapCell( :, 3) = cellfun(@(x,z)  interp2( XImgMesh, ZImgMesh, Mz0Map, x, z, 'spline'), RotMeshCell(:,3) , RotMeshCell(:,4), 'UniformOutput',false);
    % Getting the Magnetization Map(  M = Mxcos(theta) - Mzsin(theta) + iMy  )
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