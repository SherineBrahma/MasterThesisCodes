%clear all
%load('C:\Users\Sherine\Desktop\DelftStudy\Thesis\Codes\Results\AnglAlgoInvestigation.mat')
%load('C:\Users\Sherine\Desktop\DelftStudy\Thesis\Codes\Results\AngleInvestigation.mat')

LStyle = {'-','-','-','-','-',':','-.','--'};
LColor = {[0.498 0.549 0.552], [0.945 0.768 0.058], [0.152 0.682 0.376], [0.901 0.494 0.133], [0.160 0.501 0.725], [0.905 0.298 0.235], [0.556 0.266 0.678], [0.172 0.243 0.313]};
MStyle = {'.','o','s','v','x','^','d','>'};
MSize = {20,8,10,13,15,10,10,13};

AngSeqNo = 2;
NoOfInst = 3;
HblsqrSeqColtCell = {};
HtvSeqColtCell = {};
HblsqrSeqSizeMat = [];
HtvSeqSizeMat = [];
for InstCnt = 1:1:NoOfInst
    
    %IRhybrid_lsqr
    HblsqrSeqSizeCol = [];
    for AngSeqCnt = 1:1:AngSeqNo
        HblsqrSeqColtCell{InstCnt,AngSeqCnt} = StatReconsImgCell{InstCnt,1}{AngSeqCnt,2}{1,2}{1,3}.Enrm;
        Hblsqr = StatReconsImgCell{InstCnt,1}{AngSeqCnt,2}{1,2}{1,3}.Enrm;
        HblsqrSeqSizeCol = [HblsqrSeqSizeCol ; size(Hblsqr,1)];
        %figure(1), plot(Hblsqr,'color',LColor{AngSeqCnt},'LineWidth',1.5,'Marker',MStyle{AngSeqCnt},'MarkerIndices',length(Hblsqr),'MarkerSize',MSize{AngSeqCnt})
        %hold on
    end
    HblsqrSeqSizeMat = [HblsqrSeqSizeMat HblsqrSeqSizeCol];
    %hold off
    %legend('5 Angles','10 Angles','25 Angles','36 Angles','45 Angles','55 Angles','64 Angles','72 Angles')
    %legend('0-360 Deg with 10 Deg Steps','0-180 Deg with 5 Deg Steps','180-360 Deg with 5 Deg Steps','Random Selection 1','Random Selection 2','Random Selection 3','Random Selection 4','Using Coherence')
    
    %IRtv
    HtvSeqSizeCol = [];
    for AngSeqCnt = 1:1:AngSeqNo
        HtvSeqColtCell{InstCnt,AngSeqCnt} = StatReconsImgCell{InstCnt,1}{AngSeqCnt,2}{1,2}{1,5}.Enrm;
        Htv = StatReconsImgCell{InstCnt,1}{AngSeqCnt,2}{1,2}{1,5}.Enrm;
        HtvSeqSizeCol = [HtvSeqSizeCol ; size(Htv,1)];
        %figure(2), plot(Htv,'color',LColor{AngSeqCnt},'LineWidth',1.5,'Marker',MStyle{AngSeqCnt},'MarkerIndices',length(Htv),'MarkerSize',MSize{AngSeqCnt})
        %hold on
    end
    HtvSeqSizeMat = [HtvSeqSizeMat HtvSeqSizeCol];
    %hold off
    %legend('5 Angles','10 Angles','25 Angles','36 Angles','45 Angles','55 Angles','64 Angles','72 Angles')
    %legend('0-360 Deg with 10 Deg Steps','0-180 Deg with 5 Deg Steps','180-360 Deg with 5 Deg Steps','Random Selection 1','Random Selection 2','Random Selection 3','Random Selection 4','Using Coherence')
    %pause();
    
end


% Padding the sequences and truncating sequence

%IRhybrid_lsqr
HblsqrMaxSeqSize = max(HblsqrSeqSizeMat,[],2);
HblsqrMeanSeqSize = round(mean(HblsqrSeqSizeMat,2));
HblsqrPadTrunCell = {};
IntHblsqrSeqColtCell = HblsqrSeqColtCell;

%IRhtv
HtvMaxSeqSize = max(HtvSeqSizeMat,[],2);
HtvMeanSeqSize = round(mean(HtvSeqSizeMat,2));
HtvPadTrunCell = {};
IntHtvSeqColtCell = HtvSeqColtCell;
for AngSeqCnt = 1:1:AngSeqNo
    
    %IRhybrid_lsqr
    IntHblsqrPadTrunCell = cellfun(@(InstCol)  padarray(InstCol,[(HblsqrMaxSeqSize(AngSeqCnt)-size(InstCol,1)) 0],'replicate','post') , IntHblsqrSeqColtCell(:,AngSeqCnt), 'UniformOutput',false);
    IntHblsqrSeqColtCell(:,AngSeqCnt) = IntHblsqrPadTrunCell;
    IntHblsqrPadTrunCell = cellfun(@(InstCol)  TrncArry(InstCol, HblsqrMeanSeqSize(AngSeqCnt)) , IntHblsqrSeqColtCell(:,AngSeqCnt), 'UniformOutput',false);
    HblsqrPadTrunCell = [HblsqrPadTrunCell IntHblsqrPadTrunCell];
    
    %IRhtv
    IntHtvPadTrunCell = cellfun(@(InstCol)  padarray(InstCol,[(HtvMaxSeqSize(AngSeqCnt)-size(InstCol,1)) 0],'replicate','post') , IntHtvSeqColtCell(:,AngSeqCnt), 'UniformOutput',false);
    IntHtvSeqColtCell(:,AngSeqCnt) = IntHtvPadTrunCell;
    IntHtvPadTrunCell = cellfun(@(InstCol)  TrncArry(InstCol, HtvMeanSeqSize(AngSeqCnt)) , IntHtvSeqColtCell(:,AngSeqCnt), 'UniformOutput',false);
    HtvPadTrunCell = [HtvPadTrunCell IntHtvPadTrunCell];
    
end

%%%*************************************************************************************************

% Plotting the mean

%IRhybrid_lsqr
HblsqrMeanSeq = {};
HblsqrLastSeqVal = [];
for AngSeqCnt = 1:1:AngSeqNo
    HblsqrMeanSeq{1,AngSeqCnt} = mean(cell2mat(HblsqrPadTrunCell(:,AngSeqCnt)'),2);
    HblsqrMeanVal = mean(cell2mat(HblsqrPadTrunCell(:,AngSeqCnt)'),2);
    HblsqrLastSeqVal = [HblsqrLastSeqVal HblsqrMeanVal(end)];
    figure(1), plot(HblsqrMeanVal,'color',LColor{AngSeqCnt},'linestyle',LStyle{AngSeqCnt},'Marker',MStyle{AngSeqCnt},'MarkerIndices',length(HblsqrMeanVal),'MarkerSize',MSize{AngSeqCnt},'LineWidth',1.5)
    set(gcf, 'Position',  [523,420,714,389])
    %xlim([0 10])
    xlabel('No. of Iterations')
    ylabel('Relative Residual Norm')
    title('IRhybrid\_lsqr')
    hold on
end
HblsqrLastSeqVal = round(HblsqrLastSeqVal,4);
hold off
%legend('5 Angles','10 Angles','25 Angles','36 Angles','45 Angles','55 Angles','64 Angles','72 Angles')
%legend('0-360 Deg with 10 Deg Steps','0-180 Deg with 5 Deg Steps','180-360 Deg with 5 Deg Steps','Random Selection 1','Random Selection 2','Random Selection 3','Random Selection 4','Using Coherence')
legend('1','2')

%IRhtv
HtvMeanSeq = {};
HtvLastSeqVal = [];
for AngSeqCnt = 1:1:AngSeqNo
    HtvMeanSeq{1,AngSeqCnt} = mean(cell2mat(HtvPadTrunCell(:,AngSeqCnt)'),2);
    HtvMeanVal = mean(cell2mat(HtvPadTrunCell(:,AngSeqCnt)'),2);
    HtvLastSeqVal = [HtvLastSeqVal HtvMeanVal(end)];
    figure(2), plot(HtvMeanVal,'color',LColor{AngSeqCnt},'linestyle',LStyle{AngSeqCnt},'Marker',MStyle{AngSeqCnt},'MarkerIndices',length(HtvMeanVal),'MarkerSize',MSize{AngSeqCnt},'LineWidth',1.5)
    set(gcf, 'Position',  [523,420,714,389])
    xlim([0 30])
    xlabel('No. of Iterations')
    ylabel('Relative Residual Norm')
    title('IRhtv')
    hold on
end
HtvLastSeqVal = round(HtvLastSeqVal,4);
hold off
%legend('5 Angles','10 Angles','25 Angles','36 Angles','45 Angles','55 Angles','64 Angles','72 Angles')
%legend('0-360 Deg with 10 Deg Steps','0-180 Deg with 5 Deg Steps','180-360 Deg with 5 Deg Steps','Random Selection 1','Random Selection 2','Random Selection 3','Random Selection 4','Using Coherence')
legend('1','2')

