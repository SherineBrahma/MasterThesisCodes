function ReconsImg = PInvRecon(PVectMapCell, KVectMapCell, DataCell, TimeArry, StpSizeArea)

    % Get System Matrix Cell
    StoreSize = cellfun(@(MatSize)  size(MatSize) , PVectMapCell, 'UniformOutput',false);
    StkdPVectMapCell = cellfun(@(StkMat)  StkMat(:)' , PVectMapCell, 'UniformOutput',false);
    StkdKVectMapCell = cellfun(@(StkMat)  StkMat(:)' , KVectMapCell, 'UniformOutput',false);
    SysCell = cellfun(@( P, K)  StpSizeArea * (K .* exp(TimeArry' * (-conj(P)))), StkdPVectMapCell , StkdKVectMapCell, 'UniformOutput',false);
    
    %%%%%%%%%%% EXP %%%%%%%%%%%%%%%%%%

    CurrentSysangle = SysCell{2,1};
    
    AbsCurrentSysangle = abs(CurrentSysangle);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Get System Matrix
    SysMat = vertcat(SysCell{:,1});
 
    % Get Data Vector
    DataVec = vertcat(DataCell{:,1});
    
    % Apply Least Squares
    
    PInvImage = pinv(SysMat) * DataVec;
    PInvImage = abs(reshape(PInvImage, StoreSize{1,1}));
 
    % Getting the Final Image
    ReconsImg = abs(PInvImage);
    
end