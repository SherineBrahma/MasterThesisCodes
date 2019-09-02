function GradBZFit = GetGradBZFit(BZlsqCoeff, XImgCoord, ZImgCoord, StpSizeXImg, StpSizeZImg)

[XImgMesh, ZImgMesh] = meshgrid(XImgCoord, ZImgCoord);

GradBZFit = 2*abs( StpSizeXImg * BZlsqCoeff(2,1) * XImgMesh + StpSizeZImg * BZlsqCoeff(3,1) * ZImgMesh);

end