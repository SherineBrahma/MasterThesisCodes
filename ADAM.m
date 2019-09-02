function [x, ParamNormArry] = ADAM(Grad, Tolerance, InitialGuess)

    x = InitialGuess;
    m = zeros(size(InitialGuess));
    v = zeros(size(InitialGuess));
    Beta1 = .9;
    Beta2 = .999;
    eta = 0.1 ;
    epsilon = .000001;
    StpCrtraArry = zeros([10 1]);
    MaxIterations = 2000;
    DropFactor = 0.1;
    DropPeriod = 1000000;
    PrevNorm = 0;
    NoOfDigits = numel(num2str(DropPeriod));
    for ItrCount = 1:1:MaxIterations
        DscCriteria = ceil(rem(ItrCount,DropPeriod)/(10 ^ NoOfDigits));
        eta =  eta * max(DscCriteria, DropFactor - DscCriteria);
        G = - Grad(x);
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
    
end