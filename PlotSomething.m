

AngSeqNo = 4;

LStyle = {'-','-','-','-','-',':','-.','--'};
LColor = {[0.498 0.549 0.552], [0.945 0.768 0.058], [0.152 0.682 0.376], [0.901 0.494 0.133], [0.160 0.501 0.725], [0.905 0.298 0.235], [0.556 0.266 0.678], [0.172 0.243 0.313]};
MStyle = {'.','o','s','v','x','^','d','>'};
MSize = {20,8,10,13,15,10,10,13};



for i = 1:1:AngSeqNo
j = (2*(i-1))+1;
%IRhybrid_lsqr
b = StatReconsImgCell{1,1}{i,2}{1,2}{1,3}.Enrm;
plot(b,'color',LColor{j},'LineWidth',1.5,'Marker',MStyle{j},'MarkerIndices',length(b),'MarkerSize',MSize{j})
hold on
%IRtv
c = StatReconsImgCell{1,1}{i,2}{1,2}{1,5}.Enrm;
plot(c,'color',LColor{j+1},'LineWidth',1.5,'Marker',MStyle{j+1},'MarkerIndices',length(c),'MarkerSize',MSize{j+1})
hold on
end
hold off
legend('1H','1T','2H','2T','3H','3T','4H','4T')

%plot(b,'color',LColor{a},'linestyle',LStyle{a},'Marker',MStyle{a},'MarkerIndices',length(b),'MarkerSize',MSize{a},'LineWidth',1.5)