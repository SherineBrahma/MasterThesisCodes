function [ DataCell] = GetDataMatCompParams(DataStore, AnglArry)
 
    %Chopping the data sequences to include only those after the peak
    % Format of DataVect:
    % Row -> Angles Axis
    % Column -> Time Axis
    RawDataVect = (squeeze(DataStore(:,1,:)) + (1i * squeeze(DataStore(:,2,:)))).'; 
    
    %% DataCell
    % Getting DataCell
    
    DataCell = num2cell(RawDataVect, 2);
    
    DataCell = cellfun(@(x)  x.', DataCell, 'UniformOutput',false);       % Putting the same time varying data to be processed for different imaging angle rotating matrix

end