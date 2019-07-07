
function AUROC()
%%
%load data
delayLen=4;
load('FR_Trial.mat');

uidx=203;

[bA,sbA]=basestat(frA);
[bB,sbB]=basestat(frB);

delayWin=20+delayLen*10+(-29:10);
%normalization
mAZ=(mean(frA(:,delayWin),2)-bA)./sbA;
mBZ=(mean(frB(:,delayWin),2)-bB)./sbB;

%%
%plot firing rate
fh=figure('Color','w');
fh.Position(3:4)=[750,215];

subplot(1,4,1);
hold on;
cia=bootci(1000,@(x) smooth(mean(x),3),frA);
cib=bootci(1000,@(x) smooth(mean(x),3),frB);
fill([1:length(cia),fliplr(1:length(cia))],[cia(1,:),fliplr(cia(2,:))],[0.8,0.8,1],'EdgeColor','none');
fill([1:length(cib),fliplr(1:length(cib))],[cib(1,:),fliplr(cib(2,:))],[1,0.8,0.8],'EdgeColor','none');
plot(smooth(mean(frA),3),'-b','LineWidth',1.5);
plot(smooth(mean(frB),3),'-r','LineWidth',1.5);
yspan=ylim();
arrayfun(@(x) plot([x,x],yspan,':k'),[20.5,30.5,delayLen*10+30.5,delayLen*10+40.5]);
set(gca,'XTick',20.5:50:120.5,'XTickLabel',0:5:10);
xlabel('Time (s)');
ylim(yspan);
xlim([14.5,delayLen*10+41]);
ylabel(sprintf('FR (Hz) %d-%d',delayLen,uidx));

%%
%plot histogram in correct trials
subplot(1,4,2);
hold on;
bins=linspace(min([mAZ(:);mBZ(:)]),max([mAZ(:);mBZ(:)]),10);
cA=histcounts(mAZ,bins,'Normalization','probability');
cB=histcounts(mBZ,bins,'Normalization','probability');

xpos=bins(1:end-1)+diff(bins(1:2));
dpos=diff(xpos(1:2));
bar(xpos-dpos*0.225,cA,0.45,'EdgeColor','none','FaceColor','b');
bar(xpos+dpos*0.225,cB,0.45,'EdgeColor','none','FaceColor','r');
xlabel('FR Z-score');
ylabel(sprintf('Probability %dT',numel([mAZ(:);mBZ(:)])));


%%
%plot histogram in error trials
subplot(1,4,3);
hold on;
mAEZ=(mean(frAE(:,delayWin),2)-bA)./sbA;
mBEZ=(mean(frBE(:,delayWin),2)-bB)./sbB;

cAEZ=histcounts(mAEZ,bins,'Normalization','probability');
cBEZ=histcounts(mBEZ,bins,'Normalization','probability');

xpos=bins(1:end-1)+diff(bins(1:2));
dpos=diff(xpos(1:2));
bar(xpos-dpos*0.225,cAEZ,0.45,'EdgeColor','none','FaceColor','b');
bar(xpos+dpos*0.225,cBEZ,0.45,'EdgeColor','none','FaceColor','r');
xlabel('FR Z-score');
ylabel(sprintf('Probability %dT ',numel([mAEZ(:);mBEZ(:)])));


%%
%AUROC
labels=[repmat({'S1'},size(mAZ));repmat({'S2'},size(mBZ))];
scores=[mAZ;mBZ];
if mean(mAZ)>mean(mBZ)
    pos='S1';
else
    pos='S2';
end

[ax,ay,~,auc]=perfcurve(labels,scores,pos);
labelse=[repmat({'S1'},size(mAEZ));repmat({'S2'},size(mBEZ))];
scorese=[mAEZ;mBEZ];

[aex,aey,~,auce]=perfcurve(labelse,scorese,pos);
subplot(1,4,4);
hold on;
plot(ax(:,1),ay(:,1),'-b','LineWidth',1.5);text(0.33,0.66,num2str(auc(1)));
plot(aex(:,1),aey(:,1),'-r','LineWidth',1.5);text(0.66,0.33,num2str(auce(1)));
xlabel('False rate');
ylabel('True rate');

end



%%
function [mm,stdd]=basestat(fr)
frb=mean(fr(:,12:19),2);
mm=mean(frb);
stdd=std(frb);
end
