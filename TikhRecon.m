function ReconsImg = TikhRecon(Lambda, PVectMapCell, KVectMapCell, DataCell, TimeArry, StpSizeArea)

    % Get System Matrix Cell
    StoreSize = cellfun(@(MatSize)  size(MatSize) , PVectMapCell, 'UniformOutput',false);
    StkdPVectMapCell = cellfun(@(StkMat)  StkMat(:)' , PVectMapCell, 'UniformOutput',false);
    StkdKVectMapCell = cellfun(@(StkMat)  StkMat(:)' , KVectMapCell, 'UniformOutput',false);
    SysCell = cellfun(@( P, K)  StpSizeArea * (K .* exp(TimeArry' * (-conj(P)))), StkdPVectMapCell , StkdKVectMapCell, 'UniformOutput',false);
    
    % Get System Matrix
    SysMat = vertcat(SysCell{:,1});
 
    % Get Data Vector
    DataVec = vertcat(DataCell{:,1});
    
    % Get F Operator
    %F = zeros([size(SysMat,2) size(SysMat,2)]);
    %F = eye([size(SysMat,2) size(SysMat,2)]);
    %F = full(gallery('tridiag',size(SysMat,2),0,1,-1));
    F = DiffMat(StoreSize, SysMat);
    
    % Apply Tikhonov Regularization
    Lambda =  1 * 10^(-7);     % 1st 10^(-9) 2nd 10^(-7)
    ATA = (SysMat') * SysMat;
    FTF = F'*F;
    ATy =  SysMat'*DataVec;
    M = ATA + (Lambda * FTF);
    Grad = @(x) (M*x - ATy);
    TikhImage = zeros([size(SysMat,2) 1]);
    Tolerance = .000001;
    [TikhImage, ParamNormArry] = ADAM(Grad, Tolerance, TikhImage);
    TikhImage = abs(reshape(TikhImage, StoreSize{1,1}));
 
    % Getting the Final Image
    ReconsImg = abs(TikhImage);  
    
end






