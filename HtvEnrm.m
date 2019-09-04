clear all
load('C:\Users\Sherine\Desktop\DelftStudy\Thesis\Codes\Results\AngleInvestigation.mat')
%Properties -> Enrm Rnrm Xnrm
%Algorithms -> 1)IRhybrid_lsqr 2)IRhtv
AnglExpNo = 8;
LStyle = {'-','-','-','-','-',':','-.','--'};
LColor = {[0.498 0.549 0.552], [0.945 0.768 0.058], [0.152 0.682 0.376], [0.901 0.494 0.133], [0.160 0.501 0.725], [0.905 0.298 0.235], [0.556 0.266 0.678], [0.172 0.243 0.313]};
MStyle = {'.','o','s','v','x','^','d','>'};
MSize = {20,8,10,13,15,10,10,13};
A = {};
C = [];
for InsNo=1:1:20
AngoNo = 2;
B = [];
for a = 1:1:AnglExpNo
    A{InsNo,a} = StatReconsImgCell{InsNo,1}{a,2*AngoNo}.Enrm;
    IntB = size(StatReconsImgCell{InsNo,1}{a,2*AngoNo}.Enrm,1);
    B = [B ; IntB];
    %plot(StatReconsImgCell{InsNo,1}{a,2*AngoNo}.Enrm,'linestyle',LStyle{a},'LineWidth',2)
    %xlim([0 10])
    %hold on
end
C = [C B];
hold off
%legend('0-360 Deg, 10 Deg Step','0-360 Deg, 5 Deg Step','0-180 Deg, 10 Deg Step','0-180 Deg, 5 Deg Step','0-180 Deg, 2.5 Deg Step','0-90 Deg, 10 Deg Step','0-90 Deg, 5 Deg Step','0-90 Deg, 2.5 Deg Step')
%pause;
end

% Padding the sequences and truncating sequence
D = max(C,[],2);
H = round(mean(C,2));
F = {};
A1 = A;
for a = 1:1:AnglExpNo
    E = cellfun(@(InstCol)  padarray(InstCol,[(D(a)-size(InstCol,1)) 0],'replicate','post') , A1(:,a), 'UniformOutput',false);
    A1(:,a) = E;
    E = cellfun(@(InstCol)  TrncArry(InstCol, H(a)) , A1(:,a), 'UniformOutput',false);
    F = [F E];
end

% Plotting the mean
MeanSeq = {};
LastSeqVal = [];
for a = 1:1:AnglExpNo
    MeanSeq{1,a} = mean(cell2mat(F(:,a)'),2);
    b = mean(cell2mat(F(:,a)'),2);
    LastSeqVal = [LastSeqVal b(end)];
    plot(b,'color',LColor{a},'linestyle',LStyle{a},'Marker',MStyle{a},'MarkerIndices',length(b),'MarkerSize',MSize{a},'LineWidth',1.5)
    set(gcf, 'Position',  [523,420,714,389])
    xlim([0 12])
    xlabel('No. of Iterations')
    ylabel('Relative Error Norm')
    title('IRhtv')
    hold on
end
LastSeqVal = round(LastSeqVal,4);
hold off
legend('0-360 Deg, 10 Deg Step','0-360 Deg, 5 Deg Step','0-180 Deg, 10 Deg Step','0-180 Deg, 5 Deg Step','0-180 Deg, 2.5 Deg Step','0-90 Deg, 10 Deg Step','0-90 Deg, 5 Deg Step','0-90 Deg, 2.5 Deg Step')
