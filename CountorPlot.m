clear all;
x = [-10:.1:10];
y = [-10:.1:10];
[X,Y] = meshgrid(x, y);
xCoord = reshape(X,1,[]);
yCoord = reshape(Y,1,[]);
zCoord = zeros(size(yCoord));
R = rotx(20);
RotCord = R*[zCoord; xCoord; yCoord];
XRot = reshape(RotCord(2,:),size(X,1),size(X,2));
YRot = reshape(RotCord(3,:),size(Y,1),size(Y,2));
Z = sqrt(((2*XRot).^2)+((1*YRot).^2));
levels = 10:10:10;
contour( Z, levels)
axis equal
%surf(Z)
%axis equal

