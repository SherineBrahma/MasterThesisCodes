function ReconsImg = TotalVarRecon(Lambda, PVectMapCell, KVectMapCell, DataCell, TimeArry, StpSizeArea)

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
    
    % Apply Total Variation Regularization
    Lambda =  1 * 10^(-9);
    ATA = (SysMat') * SysMat;
    ATy =  SysMat'*DataVec;
    Grad = @(x) ((ATA*x - ATy) + (Lambda * F'*sign(F*x)));
    TVarImage = zeros([size(SysMat,2) 1]);
    Tolerance = .000001;
    [TVarImage, ParamNormArry] = ADAM(Grad, Tolerance, TVarImage);
    TVarImage = abs(reshape(TVarImage, StoreSize{1,1}));
 
    % Getting the Final Image
    ReconsImg = abs(TVarImage);
    figure, imagesc(ReconsImg)   
    
end