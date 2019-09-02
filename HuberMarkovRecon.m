function ReconsImg = HuberMarkovRecon(Lambda, PVectMapCell, KVectMapCell, DataCell, TimeArry, StpSizeArea)

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
    
    % Apply Huber-Markov Regularization
    Lambda1 =  1 * 10^(-9); %1 * 10^(-8);
    Lambda2 =  1 * 10^(-7); %1 * 10^(-7);
    ATA = (SysMat') * SysMat;
    ATy =  SysMat'*DataVec;
    FTF = F'*F;
    M = ATA + (Lambda1 * FTF);
    TParam = 0.005;     % Tells where the linear region start
    GradLesThnT = @(x) (M*x - ATy);
    GradGtrThnT = @(x) ((ATA*x - ATy) + (TParam * Lambda2 * F'*sign(F*x)));
    HuMkImage = zeros([size(SysMat,2) 1]);
    Tolerance = .000001;
    
    % Use ADAM Descent
    x = HuMkImage;
    m = zeros(size(x));
    v = zeros(size(x));
    G = zeros(size(x));
    Beta1 = .9;
    Beta2 = .999;
    eta = 0.1 ;
    epsilon = .000001;
    StpCrtraArry = zeros([10 1]);
    MaxIterations = 1000;
    DropFactor = 0.1;
    DropPeriod = 1000000;
    PrevNorm = 0;
    ParamNormArry = [];
    NoOfDigits = numel(num2str(DropPeriod));
    for ItrCount = 1:1:MaxIterations
        DscCriteria = ceil(rem(ItrCount,DropPeriod)/(10 ^ NoOfDigits));
        eta =  eta * max(DscCriteria, DropFactor - DscCriteria);
        GLT = -GradLesThnT(x);
        GGT = -GradGtrThnT(x);
        G(abs(F*x)<TParam) = GLT(abs(F*x)<=TParam);
        G(abs(F*x)>TParam) = GGT(abs(F*x)>TParam);
        m = (Beta1 * m) + ((1 - Beta1) * G);
        v = (Beta1 * v) + ((1 - Beta1) * (G.^ 2));
        mhat = m./(1-Beta1);
        vhat = v./(1-Beta2);
        EfStp = (eta * (mhat./((vhat).^(1/2) + epsilon)));
        x = x + EfStp;
        %ViewX = abs(reshape(x, [65 65]));
        %imagesc(ViewX)
        
        % Norm Evolution Array
        ParamNormArry(ItrCount) =  .9936 - norm(abs(x));
        plot(ParamNormArry)
        pause(0.00000001);        
        % Stopping Criteria
        StpCrtraArry = circshift(StpCrtraArry,1);
        StpCrtraArry(1) = abs(norm(abs(x)) - PrevNorm);
        PrevNorm = norm(abs(x));
        if norm(StpCrtraArry) < Tolerance
            break;
        end
    end
    
    HuMkImage = x;
    HuMkImage = abs(reshape(HuMkImage, StoreSize{1,1}));
    
    % Getting the Final Image
    ReconsImg = abs(HuMkImage);
    figure, imagesc(ReconsImg)   
    
end






