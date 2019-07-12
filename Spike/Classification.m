clear
close all
totalTrial=100;
ATrials=1:2:totalTrial;
BTrials=2:2:totalTrial;
AGrp=1:3:7;
BGrp=2:3:8;
CGrp=3:3:9;

simD=zeros(totalTrial,9);

tags=zeros(totalTrial,1);
tags(ATrials)=1;


simD(ATrials,AGrp)=rand(length(ATrials),length(AGrp)).*10+2;
simD(BTrials,AGrp)=rand(length(ATrials),length(AGrp)).*10;
simD(ATrials,BGrp)=rand(length(BTrials),length(BGrp)).*10;
simD(BTrials,BGrp)=rand(length(BTrials),length(BGrp)).*10+2;

simD(:,CGrp)=rand(totalTrial,3).*11;

figure();
imagesc(simD);
cb=colorbar();
cb.Label.String='FR (Hz)';
xlabel('Neuron #');
ylabel('Trial #');

svmModel=fitcsvm(simD(1:(totalTrial-10),:),tags(1:totalTrial-10));
classified=predict(svmModel,simD(totalTrial-9:end,:));
actual=tags((totalTrial-9):end);
fprintf('Correct Rate %d%%\n',sum(actual==classified)*10);