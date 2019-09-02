function [AcqSysMat] = GetAcqSysMat(PVectMapCell, KVectMapCell, TimeArry, StpSizeArea, AnglArry)

    % StpSizeArea area is a constant so is ingnored
    StpSizeArea = 10^(-4);
    % Get System Matrix Cell
    StoreSize = cellfun(@(MatSize)  size(MatSize) , PVectMapCell, 'UniformOutput',false);
    StkdPVectMapCell = cellfun(@(StkMat)  StkMat(:)' , PVectMapCell, 'UniformOutput',false); 
    StkdKVectMapCell = cellfun(@(StkMat)  StkMat(:)' , KVectMapCell, 'UniformOutput',false);
    %AcqSysCell = cellfun(@( P, K) StpSizeArea * (K .* (exp(TimeArry * (-P)))), StkdPVectMapCell , StkdKVectMapCell, 'UniformOutput',false);
    AcqSysCell = cellfun(@( P, K) StpSizeArea * (exp(TimeArry * (-P))), StkdPVectMapCell , StkdKVectMapCell, 'UniformOutput',false);
    % Get System Matrix
    AcqSysMat = vertcat(AcqSysCell{:,1});
    
end 