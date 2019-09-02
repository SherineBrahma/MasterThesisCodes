function [ SMat, PVectMapCell, KVectMapCell, KTpzmVectMapCell, TimeArry, StpSizeArea] = GetAcqSysMatCompParams( BZMsrd, BZMsrdIntrpol, XImgCoord, ZImgCoord, AnglArry, StepSizeTime, NoOfTimeInstances, StpSizeXImg, StpSizeZImg)

    global Gamma 
    
    %% Getting Area Step Size
    StpSizeArea = StpSizeXImg * StpSizeZImg;
    
    %% Making TimeArry and AnglArry
    % Getting the time instances used for acquisition
    TimeArry = (0 : StepSizeTime: StepSizeTime*(NoOfTimeInstances-1))';    
    % Getting the angles used for acquisition
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
    CompB0MapCell( :, 1) = cellfun(@(x,z)  interp2( XImgMesh, ZImgMesh, BZMsrdIntrpol, x, z, 'spline'), RotMeshCell(:,3) , RotMeshCell(:,4), 'UniformOutput',false);
    % Getting delta omega in the obtained rotated meshes
    BZMsrdCenter = BZMsrd(floor((size(BZMsrd,1))/2) + 1, floor((size(BZMsrd,2))/2) + 1);     % Getting the field strength of the central coordinates of the measured field
    DelOmegaMapCell = cellfun(@(B)  Gamma * (B - BZMsrdCenter) , CompB0MapCell(:,1), 'UniformOutput',false);
    %save('C:\Users\Sherine\Desktop\DelftStudy\Thesis\Codes\DataSimulationCode\DelOmegaMapCell','DelOmegaMapCell');
    %rankcell = cellfun(@(SysMatAngl)  rank(SysMatAngl) , DelOmegaMapCell, 'UniformOutput',false);
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