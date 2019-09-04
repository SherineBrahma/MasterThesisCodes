clear all;
x=-32:1:32;
y=-32:1:32;
[X,Y] = meshgrid(x,y);
Z4 = zeros(size(X));
Grad = 1.002;
B0Cntr = 0.05859;
r = [0:1:50];
for i = 1:1:(size(r,2)-1)
    Z1 = (X.^2 + Y.^2)<r(i+1)^2;
    Z2 = (X.^2 + Y.^2)<r(i)^2;
    Z3 = double(Z1-Z2);
    Z4 = Z4 + ((Grad^(i-1))*Z3);    
end
BZMsrdIntrpol = B0Cntr*Z4;
max(max(BZMsrdIntrpol))
save('ConcBMap','BZMsrdIntrpol')
imagesc(BZMsrdIntrpol)
axis equal