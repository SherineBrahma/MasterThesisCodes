function [BZFitIntrpol, BZMsrdIntrpol, XImgCoord, ZImgCoord, StpSizeXImg, StpSizeZImg] = IntrpolBZ(BZFit, BZMsrd, XMsrdCoord, ZMsrdCoord)

% Step sizes of measurement grid
StpSizeXMsrd = 50*1e-03;%*2.05;
StpSizeZMsrd = 20*1e-03;%*2.05;

% Imaging grid
RefnFactXMsrd = 32;
RefnFactZMsrd = (StpSizeZMsrd/StpSizeXMsrd) * RefnFactXMsrd;

% Step sizes of imaging grid
StpSizeXImg = StpSizeXMsrd/RefnFactXMsrd;
StpSizeZImg = StpSizeZMsrd/RefnFactZMsrd;

% Defining bounds of imaging X coordinates
XImglow = -50*1e-03;
XImghigh = 50*1e-03;

% Defining bounds of imaging Z coordinates
ZImglow = -50*1e-03;
ZImghigh = 50*1e-03;

% Imaging coordinates
XImgCoord = [XImglow: StpSizeXImg: XImghigh].';
ZImgCoord = [ZImglow: StpSizeZImg: ZImghigh].';
[XImgMesh, ZImgMesh] = meshgrid(XImgCoord, ZImgCoord);

% Interpolate to imaging grid
BZMsrdIntrpol = interp2( XMsrdCoord, ZMsrdCoord, BZMsrd, XImgMesh, ZImgMesh, 'spline');
BZFitIntrpol = interp2( XMsrdCoord, ZMsrdCoord, BZFit, XImgMesh, ZImgMesh, 'spline');

end