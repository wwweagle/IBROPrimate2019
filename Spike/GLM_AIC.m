function GLM_AIC
rpt=4;
load('GLM_Data.mat');
fillIn=@(x,y) [x(:,1:8),repmat(y,size(x,1),1)];
decay=@(x) 50-exp(-x/21.37)*50;

allTrials5=fillIn(double(allTrials5),[decay(5) 1 0 0 0 3]);
allTrials8=fillIn(double(allTrials8),[decay(8) 1 0 0 0 6]);
allTrials12=fillIn(double(allTrials12),[decay(12) 1 0 0 0 10]);
allTrialsBase=fillIn(double(allTrialsBase),[decay(5) 0 1 0 0 3]);
allTrialsNoDelay=fillIn(double(allTrialsNoDelay),[decay(0.2) 0 1 0 0 2]);
allTrialsNoLaser=fillIn(double(allTrialsNoLaser),[decay(5) 0 0 0 0 3]);
allTrials=[allTrials5;allTrials8;allTrials12;allTrialsBase;allTrialsNoDelay;allTrialsNoLaser];

matchOdor=@(x,y) ismember(x,[2 5 7])==ismember(y,[2 5 7]);

allTrials(:,end)=matchOdor(allTrials(:,1),allTrials(:,2));
allTrials(:,4)=datasample(allTrials(:,4),numel(allTrials(:,4)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% resample here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AICs=nan(rpt,7);
rsq=nan(rpt,7);
mdlAIC=cell(rpt,7);

parfor rptIdx=1:rpt
    disp(rptIdx);
    resampTrials=datasample(allTrials,length(allTrials));
    
    factors=cell(1,size(resampTrials,2));
    for i=1:length(factors)
        factors{i}=unique(resampTrials(:,i))';
    end
    
    dataMat=[];
    
    [f1,f2,f3,f4,f5,f6,f7,f10]=ndgrid(factors{1},factors{2},factors{3},factors{8},factors{9},factors{10},factors{11},factors{14});
    for idx=1:numel(f1)
        sel=all(resampTrials(:,[1:3,8:11,14])==[f1(idx),f2(idx),f3(idx),f4(idx),f5(idx),f6(idx),f7(idx),f10(idx)],2);
        if nnz(sel)>0
            perf=sum(resampTrials(sel,4))*100/nnz(sel);
            ci=bootci(rpt, @(x) sum(x)*100/length(x),resampTrials(sel,4));
            dataMat=[dataMat;perf,ci(1),ci(2),f1(idx) f2(idx) f3(idx) f4(idx) f5(idx) f6(idx) f7(idx) 0 0 f10(idx)];
        end
    end
    
    
    %     s t l g dl pd pb ps pt mt out
    %                     prevprev
    %                     corrlick
    termsMat=[0 0 0 0 0  0  0  0  0  0  0;
        1 0 0 0 0  0  0  0  0  0  0;
        0 1 0 0 0  0  0  0  0  0  0;
        %
        0 0 1 0 0  0  0  0  0  0  0;
        0 0 0 1 0  0  0  0  0  0  0;
        0 0 0 0 0  0  0  0  0  1  0;
        0 0 0 0 1  0  0  0  0  0  0;
        0 0 0 0 0  0  1  0  0  0  0;
        0 0 0 0 0  1  0  0  0  0  0;
        
        0 0 1 1 1  1  0  0  0  0  0;
        0 0 1 1 0  0  1  0  0  0  0;
        
        ];
    
    termsMatCompact=[0 0 0 0 0  0  0  0  0  0  0;
        
    0 0 0 0 0  0  0  0  0  1  0;
    0 0 0 0 1  0  0  0  0  0  0;
    0 0 1 1 1  1  0  0  0  0  0;
    ];
%     s t l g dl pd pb ps pt mt out
%                     prevprev
%                     corrlick

termsMatPBCtrl=[0 0 0 0 0  0  0  0  0  0  0;
    
0 0 0 0 0  0  0  0  0  1  0;
0 0 0 0 1  0  0  0  0  0  0;
0 0 1 1 0  0  1  0  0  0  0;
];




%     s t l g dl pd pb ps pt mt out

termsConst=[0 0 0 0 0  0  0  0  0  0  0];

termsSTM=[0 0 0 0 0  0  0  0  0  0  0;
    1 0 0 0 0  0  0  0  0  0  0;
    0 1 0 0 0  0  0  0  0  0  0;
    0 0 0 0 0  0  0  0  0  1  0];



%     s t l g dl pd pb ps pt mt out

termsSTMD=[0 0 0 0 0  0  0  0  0  0  0;
    1 0 0 0 0  0  0  0  0  0  0;
    0 1 0 0 0  0  0  0  0  0  0;
    0 0 0 0 1  0  0  0  0  0  0;
    0 0 0 0 0  0  0  0  0  1  0];


termsNoInter=[0 0 0 0 0  0  0  0  0  0  0;
    1 0 0 0 0  0  0  0  0  0  0;
    0 1 0 0 0  0  0  0  0  0  0;
    0 0 0 0 1  0  0  0  0  0  0;
    0 0 0 0 0  1  0  0  0  0  0;
    0 0 0 0 0  0  1  0  0  0  0;
    0 0 0 0 0  0  0  0  0  1  0];




termsList={termsConst,termsSTM,termsSTMD,termsNoInter,termsMat,termsMatCompact,termsMatPBCtrl};

y=double(dataMat(:,1));
X=double(dataMat(:,4:13));

for termsIdx=1:7%length(termsList)
    mdlAIC{rptIdx,termsIdx}=fitglm(X,y,termsList{termsIdx},'Categorical',[1:4,6:10],'Distribution','normal',...
        'VarNames',{'Sample','Test','Laser','Genotype','Memory_decay','Perturb_Delay','Perturb_Baseline','Perturb_Sample','Perturb_Test','Match','Correct_rate'});
    AICs(rptIdx,termsIdx)=mdlAIC{rptIdx,termsIdx}.ModelCriterion.AIC;
    rsq(rptIdx,termsIdx)=mdlAIC{rptIdx,termsIdx}.Rsquared.Ordinary;
end

end
figure('Color','w','Position',[400,400,175,210]);
hold on;
yyaxis left;

ciAIC=bootci(rpt,@(x) mean(x), AICs);
errorbar(1:7,mean(AICs),ciAIC(1,:)-mean(AICs),ciAIC(2,:)-mean(AICs),'r.','LineWidth',1);
plot(mean(AICs),'-ro','LineWidth',1);
set(gca,'YColor','k');
ylabel('Akaike information criterion','FontSize',10);
yyaxis right;

cirsq=bootci(rpt,@(x) mean(x), rsq);
errorbar(1:7,mean(rsq),cirsq(1,:)-mean(rsq),cirsq(2,:)-mean(rsq),'k.','LineWidth',1);
plot(mean(rsq),'-ko','LineWidth',1);
set(gca,'YColor','k');
ylabel('R-squared','FontSize',10);
set(gca,'XTick',1:7);
xlabel('Model #','FontSize',10,'FontName','Helvetica','Color','k');
xlim([0.5,7.5]);

savefig('AIC.fig');
print('-depsc','-painters','AIC.eps');
save('AIC.mat','AICs','rsq');

return
end