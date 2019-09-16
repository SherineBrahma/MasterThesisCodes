function GaussStdev = GaussStdev(InputVector)

InputVector = double(InputVector);
Height = hist(InputVector,size(InputVector,1))';
ZeroPadNumber = 500;
PaddedHeight = padarray(Height,[ZeroPadNumber 0],0,'both');
Incr = (max(InputVector) - min(InputVector))/size(InputVector,1);
xmin = min(InputVector) - (Incr * ZeroPadNumber);
xmax = max(InputVector) + (Incr * ZeroPadNumber);
xAxis = (xmin:Incr:xmax)';
xAxis = xAxis(1:size(PaddedHeight,1));
GaussianFit = fit(xAxis, PaddedHeight, 'gauss1');
GaussStdev = (GaussianFit.c1)/sqrt(2);

end