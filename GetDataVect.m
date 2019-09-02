function [DataVec] = GetDataVect(DataCell)

    % Get Data Matrix
    DataVec = vertcat(DataCell{:,1});


end