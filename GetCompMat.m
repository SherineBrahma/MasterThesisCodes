function [SysMat, DataVec, F] = GetCompMat(PVectMapCell, KVectMapCell, TimeArry, StpSizeArea, AnglArry, DataCell)
    
    %% ######################## SYSTEM MATRIX #######################################
    % StpSizeArea area is a constant so is ingnored
    StpSizeArea = 10^(-4);
    % Get System Matrix Cell
    StoreSize = cellfun(@(MatSize)  size(MatSize) , PVectMapCell, 'UniformOutput',false);
    StkdPVectMapCell = cellfun(@(StkMat)  StkMat(:)' , PVectMapCell, 'UniformOutput',false);
    StkdKVectMapCell = cellfun(@(StkMat)  StkMat(:)' , KVectMapCell, 'UniformOutput',false);
    %SysCell = cellfun(@( P, K)  StpSizeArea * (K .* (exp(TimeArry * (-P)))), StkdPVectMapCell , StkdKVectMapCell, 'UniformOutput',false);
    SysCell = cellfun(@( P, K)  StpSizeArea * (exp(TimeArry * (-P))), StkdPVectMapCell , StkdKVectMapCell, 'UniformOutput',false);
    
    %{
    %SysCell = cellfun(@( P, K)  StpSizeArea * (K .* exp(TimeArry * (-conj(P)))), StkdPVectMapCell , StkdKVectMapCell, 'UniformOutput',false);
    VisCell{1,1} = TimeArry;
    VisCell = repmat(VisCell, size(AnglArry,1),1);
    VisCell( :, 2) = cellfun(@( P)  abs(P), StkdPVectMapCell, 'UniformOutput',false);
    VisCell( :, 3) = cellfun(@( K)  abs(K), StkdKVectMapCell, 'UniformOutput',false);
    VisCell( :, 4) = cellfun(@( P )  TimeArry * (-P), StkdPVectMapCell, 'UniformOutput',false);
    VisCell( :, 5) = cellfun(@( P )  exp(TimeArry * (-P)), StkdPVectMapCell, 'UniformOutput',false);
    VisCell( :, 6) = cellfun(@( K )  abs(K), StkdKVectMapCell, 'UniformOutput',false);
    VisCell( :, 7) = cellfun(@( P, K )  abs(K .* exp(TimeArry * (-P))), StkdPVectMapCell, StkdKVectMapCell, 'UniformOutput',false);
    %}
    
    % Get System Matrix
    SysMat = vertcat(SysCell{:,1});
    %SysMat = zeros(size(SysCell{1,1}));
    %for i = 1:1:size(SysCell, 1)
    %    SysMat = SysMat + SysCell{i,1};
    %end
  
  %{
  %%%%%%%%%%%% EXPERIMENT %%%%%%%%%%%%%%%%%%%%%
    
  SysMat2 = vertcat(SysCell{1:18,1});
  SysMat3 = vertcat(SysCell{19:36,1});
  SysMat = SysMat2 + SysMat3;
      
  IntSumKMat = zeros(size(KVectMapCell{1,1}));
    for i = 1:1:size(KVectMapCell, 1)
        SumKMat =  IntSumKMat + KVectMapCell{i,1};
        IntSumKMat = SumKMat;
    end
  
  A = KVectMapCell{1,1} + KVectMapCell{19,1};
  
  SysMat2 = vertcat(SysCell{1:36,1});
  SysMat3 = vertcat(SysCell{37:72,1});
  SysMat = SysMat2 + SysMat3;
  
  % Visualizing changes in the stacked magneticfiueld
  A =   SysCell{1,1};
  AStkd = StkdPVectMapCell{1,1};
  RankA = rank(A);
  RankSysMat = rank(SysMat);
  RankSysMat1 = rank(SysMat1);
  
  for i =1:1:36
  AStkd2 = StkdPVectMapCell{i,1};
  imagesc(abs(AStkd2)')
  pause(.1)
  end
  
  % Investigating independency of the columns
  
  A = vertcat(SysCell{1:3,1});
  RankA = rank(A);
  [R,basiccol] = rref(A);
  corrA = corrcoef(A);
  abscorrA = abs(corrA);
  corrAT = corrcoef(A');
  abscorrAT = abs(corrA');
  [M, MPos] = max(abscorrA,[],1);
  
  Exp = sort(abscorrAT(abscorrAT<1));
  [U,S,V] = svd(A);
  D = diag(S);
  		
  % Adding Rows to remove angle dependencies
  
  SysMat2 = vertcat(SysCell{1:36,1});
  SysMat3 = vertcat(SysCell{37:72,1});
  SysMat = SysMat2 + SysMat3;
  
  %39 column
  a = exp(-(-1.11111111111111e-16 + 10.8595865838772i));
  
  %50 column
  b = exp(-(-1.11111111111111e-16 + 4.57643195755026i));
  
  
  X = [1 2;0 0; 1 9; 0 0; 1 2];
  corrX = corrcoef(X);
  
  
  
  
  
  B = A(:,basiccol);
    
  UnqCol = unique(AStkd);
  
  B = sum(A,2);
  
  A = abs(A);
  
  
  
  
  A = PVectMapCell{:,1};
  b = abs(A);
  imagesc(b)
    
  [a, b] = min(min(b))

    CondNo4 = cond(SysMat4);
    Rank4 = rank(SysMat4);   
    
    
    %CondNo = cond(SysMat);
    %CondNo1 = cond(SysMat1);
    
    Rank = rank(SysMat);
    Rank1 = rank(SysMat1);
    %}
    
    % Rank for individual matrices for corresponding angles and Matrix used
    % for computations
    %rankcell = cellfun(@(SysMatAngl)  rank(SysMatAngl) , SysCell, 'UniformOutput',false);
    %RankSysMat = rank(SysMat);
    
    %SysMat = SysMat1;
    
    %save("4000TimeSamp", "rankcell", "RankSysMat");
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