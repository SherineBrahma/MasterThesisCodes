function ReconsImg = AdjOpRecon(PVectMapCell, KTpzmVectMapCell, DataCell, SMat, TimeArry, AnglArry, StepSizeTime, StepSizeAngl)

    % Sum them along time axis
    TimeSummedImgCell = cellfun(@(B)  zeros(size(B)) , PVectMapCell, 'UniformOutput',false);
    for TimeIndx = 1: 1: size(TimeArry,2)
            
            IntTimeSummedImgCell = cellfun(@( P, K, D)  StepSizeAngl * StepSizeTime * (K .* exp(-conj(P)*TimeArry(TimeIndx))) * SMat(TimeIndx, TimeIndx) * D(TimeIndx), PVectMapCell , KTpzmVectMapCell, DataCell, 'UniformOutput',false);
            TimeSummedImgCell = cellfun(@( TSI, InTSI)  TSI + InTSI, TimeSummedImgCell, IntTimeSummedImgCell, 'UniformOutput',false);
            %imagesc(abs(TimeSummedImgCell{1,1}))
            %pause;
            
    end
    
    % Sum them along angle axis
    AnglSummedImg = zeros(size(PVectMapCell{1,1}));
    for AnglIndx = 1: 1: size(AnglArry,1)
        
        IntAnglSummedImg = TimeSummedImgCell{ AnglIndx, 1};
        AnglSummedImg = AnglSummedImg + IntAnglSummedImg;
        
    end
 
    % Getting the Final Image
    ReconsImg = abs(AnglSummedImg);
    
end