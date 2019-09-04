clear all
load('C:\Users\Sherine\Desktop\DelftStudy\Thesis\Codes\Results\AngleInvestigation.mat')
%Properties -> Time Elapsed
%Algorithms -> 1)IRhybrid_lsqr 2)IRhtv
AnglExpNo = 8;
BColor = {'b','g','r','c','m','y','k','w'};
B = [];
for InsNo=1:1:20
AngoNo = 2;
A = [];
for a = 1:1:AnglExpNo
    IntA = StatReconsImgCell{InsNo,1}{a,2*AngoNo}.TimeElaspsed;
    A = [A IntA];
    %bar([a], StatReconsImgCell{InsNo,1}{a,2*AngoNo}.TimeElaspsed,BColor{a})
    hold on
end
B = [B; A];
hold off
%legend('0-360 Deg, 10 Deg Step','0-360 Deg, 5 Deg Step','0-180 Deg, 10 Deg Step','0-180 Deg, 5 Deg Step','0-180 Deg, 2.5 Deg Step','0-90 Deg, 10 Deg Step','0-90 Deg, 5 Deg Step','0-90 Deg, 2.5 Deg Step')
%pause;
end
C = mean(B, 1);
D = round(C,4);
for a = 1:1:AnglExpNo
    bar([a],C(a), BColor{a})
    hold on
end
hold off
legend('0-360 Deg, 10 Deg Step','0-360 Deg, 5 Deg Step','0-180 Deg, 10 Deg Step','0-180 Deg, 5 Deg Step','0-180 Deg, 2.5 Deg Step','0-90 Deg, 10 Deg Step','0-90 Deg, 5 Deg Step','0-90 Deg, 2.5 Deg Step')
