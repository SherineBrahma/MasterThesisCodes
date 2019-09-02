function [BZMsrd, XMsrdCoord, ZMsrdCoord] = LoadMsrdBZField(MsrdBZDataFile)

BMsrdParams = load(MsrdBZDataFile);

BZMsrd = BMsrdParams.Bz;           % Measured Z component of Magnetic Field, (ZCoord, XCoord)
XMsrdCoord = BMsrdParams.Xc;       % X coordinates of the measurement grid
ZMsrdCoord = BMsrdParams.Zc;       % Z coordinates of the measurement grid

end