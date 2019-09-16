clear all;
x=-32:1:32;
y=-32:1:32;
[X,Y] = meshgrid(x,y);
B0Cntr = 0.05859;
Grad = .0083;
Z4 = B0Cntr +  Grad * (X./65);
BZMsrdIntrpol = Z4;
%save('GradBMap','BZMsrdIntrpol')
imagesc(BZMsrdIntrpol)
axis equal