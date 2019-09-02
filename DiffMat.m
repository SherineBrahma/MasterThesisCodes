function DiffMat = DiffMat(StoreSize, SysMat)

    m = StoreSize{1,1}(1);
    n = StoreSize{1,1}(2);
    ZeroMat = zeros([m + 2, n + 2]);
    %DiffKernel = [1 -2 1; -2 4 -2;1 -2 1];
    DiffKernel = [-1 -1 -1; -1 8 -1;-1 -1 -1];
    DiffMat = zeros([size(SysMat,2) size(SysMat,2)]);
    for RowIndx = 1:1:m
        for ColIndx = 1:1:n
            IntDiffMat = ZeroMat;
            
            IntDiffMat(RowIndx + 0, ColIndx + 0) = DiffKernel(1,1);
            IntDiffMat(RowIndx + 0, ColIndx + 1) = DiffKernel(1,2);
            IntDiffMat(RowIndx + 0, ColIndx + 2) = DiffKernel(1,3);
            
            IntDiffMat(RowIndx + 1, ColIndx + 0) = DiffKernel(2,1);
            IntDiffMat(RowIndx + 1, ColIndx + 1) = DiffKernel(2,2);
            IntDiffMat(RowIndx + 1, ColIndx + 2) = DiffKernel(2,3);
            
            IntDiffMat(RowIndx + 2, ColIndx + 0) = DiffKernel(3,1);
            IntDiffMat(RowIndx + 2, ColIndx + 1) = DiffKernel(3,2);
            IntDiffMat(RowIndx + 2, ColIndx + 2) = DiffKernel(3,3);
            
            IntDiffMat(1,:) = [];
            IntDiffMat(size(IntDiffMat,1),:) = [];
            IntDiffMat(:,1) = [];
            IntDiffMat(:,size(IntDiffMat,2)) = [];
            
            DiffMat((((RowIndx - 1) * n ) + ColIndx),:) = IntDiffMat(:)'; 
            
        end
    end     
end