function Selectivity
close all;
load('FR_Trial.mat');
MeanSelect=(mean(frB)-mean(frA))./(mean(frB)+mean(frA));
figure('Color','w','Position',[100,100,400,400]);
hold on;
plot((-2:0.1:10.9)+0.05,smooth(MeanSelect),'-k','LineWidth',2);
arrayfun(@(x) plot([x,x],ylim(),'--k'),[0,1,5,6]);
xlim([-1,7]);
ylim([-0.5,1]);
xlabel('Time (s)');
ylabel('Selectivity');

end