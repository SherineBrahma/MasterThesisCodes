function M0MapStkdNoisy = AddNoiseM0Map( M0MapStkd, AngleNoise)

    M0MapStkdNoisy = [];
    for PixCnt = 1:1:size(M0MapStkd,2)
        LowLim = AngleNoise(1);
        TopLim = AngleNoise(2);
        DelAlpX = LowLim + (TopLim - LowLim)*rand(1,1);
        DelAlpY = LowLim + (TopLim - LowLim)*rand(1,1);
        DelAlpZ = LowLim + (TopLim - LowLim)*rand(1,1);
        Rx = [1 0 0;0 cosd(DelAlpX) sind(DelAlpX); 0 -sind(DelAlpX) cosd(DelAlpX)];
        Ry = [cosd(DelAlpY) 0 -sind(DelAlpY);0 1 0; sind(DelAlpY) 0 cosd(DelAlpY)];
        Rz = [cosd(DelAlpX) sind(DelAlpZ) 0;-sind(DelAlpZ) cosd(DelAlpZ) 0; 0 0 1];
        M0MapStkdNoisy = [M0MapStkdNoisy Rz*Ry*Rx*M0MapStkd(:,PixCnt)];
    end
    
end