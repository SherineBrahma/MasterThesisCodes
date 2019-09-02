function [BZFit, BZlsqCoeff] = QuadFitBZMsrd(BZMsrd, XMsrdCoord, ZMsrdCoord)

BZMsrdTrnspd = BZMsrd.';           % Transposing BZMsrd, (XCoord, ZCoord)
bZMsrd = BZMsrdTrnspd(:);           % Stacking BZMsrdTrnspd to 1D array

% Trying to get in a form of -
% y = [1 + x^2 + y^2] [a b c].' => y = Blsq Coeff.'
BZlsq = zeros(size(bZMsrd,1),3);                    % Defining a Least Square System Matrix
BZlsq(:,1) = ones(size(bZMsrd,1),1);                % First column
BZlsq(:,2) = kron(ones(size(BZMsrd,1),1), XMsrdCoord.^2);   % Second column, Stacking is done in X first
BZlsq(:,3) = kron( ZMsrdCoord.^2,ones(size(BZMsrd,2), 1));   % Third column

BZlsqCoeff = ((BZlsq.' * BZlsq)^(-1)) * BZlsq.' * bZMsrd; % Getting least-squares coeffcients

[XMsrdMesh, ZMsrdMesh] = meshgrid( XMsrdCoord, ZMsrdCoord);  % Getting a Meshgrid of the XMsrdCoord and ZMsrdCoord

BZFit = BZlsqCoeff(1,1) + BZlsqCoeff(2,1) * XMsrdMesh.^2 + BZlsqCoeff(3,1) * ZMsrdMesh.^2; % Quadratic Profile

end