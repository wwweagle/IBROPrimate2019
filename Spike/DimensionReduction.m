close all;
clear;
%%
%generate simulated data
%ramping component, 8/14
a1=0.1:0.1:10;
a2=a1.*2;
a3=a1.*3;
a4=a1+1;
a5=a1+2;
a6=2*a1+10;
a7=-2*a1+20;
a8=-4*a1+40;
%wave component
s1=10*sin(linspace(0,8*pi,100))+10;
s2=0.5*s1+10;
s3=-0.5*s1+10;
s4=-s1+20;
%random noise
r1=rand(1,100)*5;
r2=rand(1,100)*5;

%combined
tmat=[a1',a2',a3',a4',a5',a6',a7',a8',s1',s2',s3',s4',r1',r2'];
%standard normalization
tmatNorm=(tmat-mean(tmat))./std(tmat);

%%
%component analysis
[coeff,score,latent]=pca(tmat);

%visualization
fh=figure();
fh.WindowState='maximized';
drawnow();
subplot(2,6,1);
plot(tmat);
xlabel('Time (s)');
ylabel('Firing rate (Hz)');
for idx=1:4
    subplot(2,6,1+idx);
    plot(score(:,idx));
    text(10,max(ylim())-0.1*(diff(ylim())),sprintf('PC%d, %2.0f%% Var',idx,latent(idx)*100/sum(latent)));
end
subplot(2,6,6);
imagesc(coeff(:,1:4),[-0.5,0.5]);
colormap('jet');colorbar();

%repeat on normalized data
[coeff,score,latent]=pca(tmatNorm);
subplot(2,6,7);
plot(tmatNorm);
xlabel('Time (s)');
ylabel('Normalized FR');
for idx=1:4
    subplot(2,6,7+idx);
    plot(score(:,idx));
    text(10,max(ylim())-0.1*(diff(ylim())),sprintf('PC%d, %2.0f%% Var',idx,latent(idx)*100/sum(latent)));
end

subplot(2,6,12);
imagesc(coeff(:,1:4),[-0.5,0.5]);
colormap('jet');colorbar();
