function BehavioralRegression
fstr=load('VarLen.mat');
dataset={'out5','out8','out12','out16','out20','out30','out40'};
measure={'perf','false','miss','dpc','lickEff'};
mDesc={'Correct rate','False choice','Miss','D prime','Lick efficiency'};
lenX=[5,8,12,16,20,30,40];

%%%%%%%%%%%COLOR MAP%%%%%%%%%%%%%%
cmap=colormap('jet');
subCmap=cmap(ceil(rand(19,1)*length(cmap)),:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

miceId=unique(cell2mat(cellfun(@(x) fstr.(x).perf(:,1),dataset,'UniformOutput',false).'));

close all;
allData=cell(1,7);
for mIdx=1%:length(measure)
    m=measure{mIdx};
    for len=1:length(dataset)
        data=fstr.(dataset{len}).(m);
        allData{len}=data(data(:,2)==0,[3 1]);
    end
    ft = fittype('50+b*exp(a*x)');
    fitData=cell2mat(arrayfun(@(x) [repmat(lenX(x),length(allData{x}),1),allData{x}(:,1)],1:length(lenX),'UniformOutput',false)');
    [fitM,gof]=fit(fitData(:,1),(fitData(:,2)),ft,'StartPoint',[-0.5,50],'Upper',[0,50]);
    
    figure('Color','w','Position',[mIdx*400+100,100,350,240]);
    hold on;
   
    for i=1:7
        arrayfun(@(x) plot(lenX(i)+rand()*1.5-0.75,allData{i}(x,1),'o','MarkerSize',4,'MarkerFaceColor',subCmap(miceId==allData{i}(x,2),:),'MarkerEdgeColor','none'),1:length(allData{i}));
        fprintf('%d s, n = %d\n',lenX(i),length(allData{i}));
    end
    plot(lenX,exp(fitM.a*lenX)*50+50,'--k','LineWidth',2);
    text(10,100,sprintf('\\tau = %0.2f (s)',-1/fitM.a),'FontSize',12,'Interpreter','tex');
    text(20,90,sprintf('r^2 = %0.3f',gof.rsquare),'FontSize',12);
    set(gca,'YTick',[50 75 100],'FontSize',12);
    xlim([0,41]);
    ylabel(mDesc{mIdx},'FontSize',12,'FontName','Helvetica');
    xlabel('Delay duration (s)','FontSize',12,'FontName','Helvetica');
end

end