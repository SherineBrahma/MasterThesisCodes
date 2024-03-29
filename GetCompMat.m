function [SysMat, DataVec, F] = GetCompMat(PVectMapCell, KVectMapCell, TimeArry, StpSizeArea, AnglArry, DataCell)
    
    %% ######################## SYSTEM MATRIX #######################################
    % StpSizeArea area is a constant so is ingnored
    StpSizeArea = 10^(-4);
    % Get System Matrix Cell
    StoreSize = cellfun(@(MatSize)  size(MatSize) , PVectMapCell, 'UniformOutput',false);
    StkdPVectMapCell = cellfun(@(StkMat)  StkMat(:).' , PVectMapCell, 'UniformOutput',false);
    StkdKVectMapCell = cellfun(@(StkMat)  StkMat(:).' , KVectMapCell, 'UniformOutput',false);
    SysCell = cellfun(@( P, K)  StpSizeArea * (K .* (exp(TimeArry * (-P)))), StkdPVectMapCell , StkdKVectMapCell, 'UniformOutput',false);
    %SysCell = cellfun(@( P, K)  StpSizeArea * (-1i .* (exp(TimeArry * (-P)))), StkdPVectMapCell , StkdKVectMapCell, 'UniformOutput',false);
    
    % Get System Matrix
    SysMat = vertcat(SysCell{:,1});
    %% ######################## DATA MATRIX #######################################
    
    % Get Data Matrix
    DataVec = vertcat(DataCell{:,1});
    %DataVec = zeros(size(DataCell{1,1}));
    %for i = 1:1:size(DataCell, 1)
    %    DataVec = DataVec + DataCell{i,1};
    %end
    
    %% ######################## REGULARIZATION MATRIX #######################################
    % Get Regularization Matrix
    F = DiffMat(StoreSize, SysMat);

end