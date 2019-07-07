function PSTH
close all;
clear;
fn='Pir51-4415-4s-day3.mat';
delayLen=4;
tetidx=8;
uidx=1;

binSize=0.25;
javaaddpath('spk2fr.jar');
s2f=spk2fr.Spk2fr;
s2f.setRefracRatio(0.1);
s2f.setLeastFR('all');
rLim=delayLen+2;% 9 for 4s delay, 13 for 8s delay

ft=load(fn);
spk=ft.Spk;

figure('Color','w','Position',[100,100,400,400]);
mm=mean(spk(all(spk(:,1:2)==[tetidx,uidx],2),4:end));
mm((1:4)*52)=nan;
plot(mm,'k-','LineWidth',2);
xlim([0,208]);

info=ft.TrialInfo;
genOne(spk,info,true(size(info,1),1),true,tetidx,uidx);


    function genOne(spk,info,filter,isCorrect,tetidx,uidx)
        ts=s2f.getTS(info(filter,:),spk,'wjdnms',false,isCorrect);
        keys=s2f.getKeyIdx();
        for k=1:size(keys,1)
            if isequal([tetidx,uidx],keys(k,:)) && (~isempty(k))
                plotOne(ts{k});
                assignin('base','TimeStamp',ts{k});
            end
        end
    end


    function plotOne(ts)
        figure('Color','w','Position',[100,100,400,400]);
        
        subplot('Position',[0.1,0.5,0.85,0.40]);
        hold on;
        
        pos=round((length(ts{1})-20)/2);
        posRange=pos:pos+19;
        cellfun(@(x) plot([x{1}';x{1}'],repmat([x{2};x{2}+0.8]+1,1,length(x{1})),'-r'),cellfun(@(x,y) {x,y},ts{1}(posRange),num2cell(1:length(posRange),1)','UniformOutput',false));
        
        pos=round((length(ts{2})-20)/2);
        posRange=pos:pos+19;
        cellfun(@(x) plot([x{1}';x{1}'],repmat([x{2};x{2}+0.8]+1,1,length(x{1})),'-b'),cellfun(@(x,y) {x,y},ts{2}(posRange),num2cell([1:length(posRange)]+20,1)','UniformOutput',false));
        
        xlim([-1,rLim]);
        ylim([0,41]);
        set(gca,'XTick',[],'YTick',[0,18,23,40],'YTickLabel',[0,20,0,20]);
        arrayfun(@(x) plot([x,x],ylim(),'--k'),[0,1,1,2,3,3.5]+[0,0,ones(1,4).*delayLen]);
        ylabel('Trial #');
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        subplot('Position',[0.1,0.1,0.85,0.35]);
        hold on;
        
        pfHist=cell2mat(cellfun(@(x) histcounts(x,-4:binSize:rLim)./binSize,ts{1},'UniformOutput',false));
        cia=bootci(1000,@(x) mean(x), pfHist);
        fill([-4+binSize/2:binSize:rLim,rLim-binSize/2:-binSize:-4],[cia(1,:),fliplr(cia(2,:))],[1,0.8,0.8],'EdgeColor','none');
        
        bnHist=cell2mat(cellfun(@(x) histcounts(x,-4:binSize:rLim)./binSize,ts{2},'UniformOutput',false));
        cib=bootci(1000,@(x) mean(x), bnHist);
        fill([-4+binSize/2:binSize:rLim,rLim-binSize/2:-binSize:-4],[cib(1,:),fliplr(cib(2,:))],[0.8,0.8,1],'EdgeColor','none');
        
        
        plot(-4+binSize/2:binSize:rLim,(mean(pfHist,1))','-r');
        plot(-4+binSize/2:binSize:rLim,(mean(bnHist,1))','-b');
        
        xlim([-1,rLim]);
        ylim([min(ylim()),max([cia(:);cib(:)])]);
        
        arrayfun(@(x) plot([x,x],ylim(),'--k'),[0,1,1,2,3,3.5]+[0,0,ones(1,4).*delayLen]);
        set(gca,'XTick',[0,5,10]);
        
        xlabel('Time(s)');
        ylabel('FR (Hz)');
    end
end



