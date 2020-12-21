
%% Analysis sketches
clear all      
cd('T:\jan\Collabo Data\HPCpaperPreProcessed')
load AATC_Sua_Psth_1ms
% load SuaGlMpValues.mat

idx=find(diff(CellPerSesCounter)~=1)
mean(CellPerSesCounter(idx))

cd('T:\jan\Collabo Data\PFCpaperPreProcessed')
load('LickEvokedIndx.mat')
%% Trace Down
Bin=25;
window=3000/Bin:4000/Bin;
baseline=1:1000/Bin;
time=-1+0.001*Bin:.001*Bin:4;

AATC_Sua_PsthBined=squeeze(mean(reshape(AATC_Sua_Psth,Bin,size(AATC_Sua_Psth,1)/Bin,size(AATC_Sua_Psth,2),size(AATC_Sua_Psth,3))));
%Subtract Baseline and Normalizse

bc_Psths=(AATC_Sua_PsthBined-nanmean(AATC_Sua_PsthBined(baseline,:,:)));%./nanmean(AATC_Sua_PsthBined(baseline,:,:))*100;;  %in herz


for i=1:size(bc_Psths,2)
   bc_Psths(:,i,:) = smoothdata(bc_Psths(:,i,:),'gaussian',25);
end

Tresh=-1;
BaseStd=squeeze(std(bc_Psths(baseline,:,:)));
EvR=squeeze(mean(bc_Psths(window,:,:)));
EvokedPeaks=squeeze(mean(bc_Psths(window,:,:)));
Evokedup(1:size(AATC_Sua_Psth,2))=0;
Evokedup(find(EvR(:,1)<BaseStd(:,1)*Tresh|EvR(:,2)<BaseStd(:,2)*Tresh))=1;

TraceDownHPC=Evokedup;

figure()
plotwindow=1:5000/Bin;
Condition1=find(LearnedCounter==1&Evokedup==1&LickUpHPC==0&LickDownHPC==0);
Condition=Condition1 

subplot(2,3,[1,2,4,5])

hold on
pp(1,:)=nanmean(bc_Psths(:,Condition,1)');
errBar=repmat(std(bc_Psths(:,Condition,2)')/sqrt(size(Condition,2)),2,1);
shadePlot2(time(plotwindow),pp(1,plotwindow),errBar(:,plotwindow),[0 0 1],.3)
hold on
pp(2,:)=nanmean(bc_Psths(:,Condition,2)');
errBar=repmat(std(bc_Psths(:,Condition,2)')/sqrt(size(Condition,2)),2,1);
shadePlot2(time(plotwindow),pp(2,plotwindow),errBar(:,plotwindow),[1 0 0],.3)
hold on
plot([0 0],[max(max(pp))+max(max(errBar)) min(min(pp))-min(min(errBar))],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
plot([2 2],[max(max(pp))+max(max(errBar)) min(min(pp))-min(min(errBar))],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
plot([3 3],[max(max(pp))+max(max(errBar)) min(min(pp))-min(min(errBar))],'LineWidth',2,'LineStyle',':','Color',[0 0 0])
box off
axis tight
xlabel('Time [s]')
ylabel('{\Delta} FiringRate [hz]')
set(gca,'FontSize',25);
set(gca,'LineWidth',5);


subplot(2,3,[3,6])
data{1}=EvokedPeaks(Condition,1);
data{2}=EvokedPeaks(Condition,2);
[p1,h]=ranksum(data{1},data{2})
y=[mean(data{1}),mean(data{2})]
e=[std(data{1})/sqrt(length(data{1})),...
std(data{2})/sqrt(length(data{2}))];
p=[ NaN p1;
    p1,NaN]
c=[0 0 1;1 0 0]
superbar(y,'E',e,'P',p,'BarFaceColor',c)
names={'CS+','CS-'}
set(gca,'Xtick',0:1:3)
set(gca,'xticklabel',names)
legend boxoff

% subplot(2,4,[8])
% NT(1)=length(find(Evokedup==1&neuronType==1&LearnedCounter==1))/length(find(neuronType==1&LearnedCounter==1))
% NT(2)=length(find(Evokedup==1&neuronType==2&LearnedCounter==1))/length(find(neuronType==2&LearnedCounter==1))
% bar(NT)
% 
% subplot(2,4,[8])
% NT(1)=length(find(Evokedup==1&neuronType==1&LearnedCounter==1))/length(find(Evokedup==1&LearnedCounter==1))
% NT(2)=length(find(Evokedup==1&neuronType==2&LearnedCounter==1))/length(find(Evokedup==1&LearnedCounter==1))
% NT(3)=length(find(Evokedup==1&(neuronType==0|neuronTypeGmm==3)&LearnedCounter==1))/length(find(Evokedup==1&LearnedCounter==1))
% pie(NT,[1 2 3])



clearvars data
Condition1=find(LearnedCounter==1&Evokedup==1&LickUpHPC==0&LickDownHPC==0);
Condition2=find(LearnedCounter==0&TrgDayCounter<3&Evokedup==1);

data{1}=EvokedPeaks(Condition1,1)';
data{2}=EvokedPeaks(Condition1,2)';
data{3}=EvokedPeaks(Condition2,1)';
data{4}=EvokedPeaks(Condition2,2)';

Reward=zeros(size([data{1},data{2},data{3},data{4}],2),1);
Reward([1:size(data{1},2),(size(data{1},2)*2)+1:size(data{1},2)*2+1+size(data{3},2)])=1;
Learned=zeros(size([data{1},data{2},data{3},data{4}],2),1);
Learned(1:size(data{1},2)+size(data{2},2))=1;
% [h,p]=anovan([data{1},data{2},data{3},data{4}],{Reward;Learned},'model','interaction')


AllPostCells=length(find(LearnedCounter==1))
PostTraceDown=length(find(LearnedCounter==1&Evokedup==1&LickUpHPC==0&LickDownHPC==0))
PerPostTracedown=PostTraceDown/AllPostCells
% figure()
% scatter(LickRateCounter(Evokedup==1),EvokedPeaks(Evokedup==1,1)-EvokedPeaks(Evokedup==1,2))
% lsline
% axis tight
% 
% [h,p]=corrcoef(LickRateCounter(Evokedup==1),EvokedPeaks(Evokedup==1,1)-EvokedPeaks(Evokedup==1,2))
% 
% load SuaGLMpValues.mat
% LickMod=length(find(SuaGLMpValues(Condition1,8)<0.05)==1)/length(Condition1)

%% Trace Up
Bin=25;
window=3000/Bin:4000/Bin;
baseline=1:1000/Bin;
time=-1+0.001*Bin:.001*Bin:4;

AATC_Sua_PsthBined=squeeze(mean(reshape(AATC_Sua_Psth,Bin,size(AATC_Sua_Psth,1)/Bin,size(AATC_Sua_Psth,2),size(AATC_Sua_Psth,3))));
%Subtract Baseline and Normalizse

bc_Psths=((AATC_Sua_PsthBined)-nanmean(AATC_Sua_PsthBined(baseline,:,:)));%./nanmean(AATC_Sua_PsthBined(baseline,:,:))*100;  %in herz


for i=1:size(bc_Psths,2)
   bc_Psths(:,i,:) = smoothdata(bc_Psths(:,i,:),'gaussian',25);
end


Tresh=1;
BaseStd=squeeze(std(bc_Psths(baseline,:,:)));
EvR=squeeze(mean((bc_Psths(window,:,:))));
EvokedPeaks=squeeze(mean(bc_Psths(window,:,:)));
Evokedup(1:size(AATC_Sua_Psth,2))=0;
Evokedup(find(EvR(:,1)>BaseStd(:,1)*Tresh|EvR(:,2)>BaseStd(:,2)*Tresh))=1;
MeanBaselineX=squeeze(mean(AATC_Sua_PsthBined(baseline,:,1)));

TraceUpHPC=Evokedup;

figure()

plotwindow=1:5000/Bin;
Condition1=find(LearnedCounter==1&Evokedup==1&LickUpHPC==0&LickDownHPC==0);
Condition=Condition1
mean(FRperiodsCounter(LearnedCounter==1&Evokedup==0,1))
subplot(2,3,[1,2,4,5])

hold on
pp(1,:)=nanmean(bc_Psths(:,Condition,1)');
errBar=repmat(std(bc_Psths(:,Condition,2)')/sqrt(size(Condition,2)),2,1);
shadePlot2(time(plotwindow),pp(1,plotwindow),errBar(:,plotwindow),[0 0 1],.3)
hold on
pp(2,:)=nanmean(bc_Psths(:,Condition,2)');
errBar=repmat(std(bc_Psths(:,Condition,2)')/sqrt(size(Condition,2)),2,1);
shadePlot2(time(plotwindow),pp(2,plotwindow),errBar(:,plotwindow),[1 0 0],.3)
hold on
plot([0 0],[max(max(pp))+max(max(errBar)) min(min(pp))-min(min(errBar))],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
plot([2 2],[max(max(pp))+max(max(errBar)) min(min(pp))-min(min(errBar))],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
plot([3 3],[max(max(pp))+max(max(errBar)) min(min(pp))-min(min(errBar))],'LineWidth',2,'LineStyle',':','Color',[0 0 0])
box off
axis tight
xlabel('Time [s]')
ylabel('{\Delta} FiringRate [hz]')
set(gca,'FontSize',25);
set(gca,'LineWidth',5);


subplot(2,3,[3,6])
data{1}=EvokedPeaks(Condition,1);
data{2}=EvokedPeaks(Condition,2);
[p1,h]=ranksum(data{1},data{2})
y=[mean(data{1}),mean(data{2})]
e=[std(data{1})/sqrt(length(data{1})),...
std(data{2})/sqrt(length(data{2}))];
p=[ NaN p1;
    p1,NaN]
c=[0 0 1;1 0 0]
superbar(y,'E',e,'P',p,'BarFaceColor',c)
names={'CS+','CS-'}
set(gca,'Xtick',0:1:3)
set(gca,'xticklabel',names)
legend boxoff


% subplot(2,4,[8])
% NT(1)=length(find(Evokedup==1&neuronType==1&LearnedCounter==1))/length(find(neuronType==1&LearnedCounter==1))
% NT(2)=length(find(Evokedup==1&neuronType==2&LearnedCounter==1))/length(find(neuronType==2&LearnedCounter==1))
% bar(NT)
% subplot(2,4,[8])
% NT(1)=length(find(Evokedup==1&neuronType==1&LearnedCounter==1))/length(find(Evokedup==1&LearnedCounter==1))
% NT(2)=length(find(Evokedup==1&neuronType==2&LearnedCounter==1))/length(find(Evokedup==1&LearnedCounter==1))
% NT(3)=length(find(Evokedup==1&(neuronType==0|neuronTypeGmm==3)&LearnedCounter==1))/length(find(Evokedup==1&LearnedCounter==1))
% pie(NT,[1 2 3])



clearvars data
Condition1=find(LearnedCounter==1&Evokedup==1&LickUpHPC==0&LickDownHPC==0);
Condition2=find(LearnedCounter==0&TrgDayCounter<3&Evokedup==1);

data{1}=EvokedPeaks(Condition1,1)';
data{2}=EvokedPeaks(Condition1,2)';
data{3}=EvokedPeaks(Condition2,1)';
data{4}=EvokedPeaks(Condition2,2)';

Reward=zeros(size([data{1},data{2},data{3},data{4}],2),1);
Reward([1:size(data{1},2),(size(data{1},2)*2)+1:size(data{1},2)*2+1+size(data{3},2)])=1;
Learned=zeros(size([data{1},data{2},data{3},data{4}],2),1);
Learned(1:size(data{1},2)+size(data{2},2))=1;
% [h,p]=anovan([data{1},data{2},data{3},data{4}],{Reward;Learned},'model','interaction')

AllPostCells=length(find(LearnedCounter==1))
PostTraceUp=length(find(LearnedCounter==1&Evokedup==1&LickUpHPC==0&LickDownHPC==0))
PerPostTraceUp=PostTraceUp/AllPostCells
% figure()
% scatter(LickRateCounter(Evokedup==1),EvokedPeaks(Evokedup==1,1)-EvokedPeaks(Evokedup==1,2))
% lsline
% axis tight
% 
% [h,p]=corrcoef(LickRateCounter(Evokedup==1),EvokedPeaks(Evokedup==1,1)-EvokedPeaks(Evokedup==1,2))
% 
% 
% %%
% 
% load SuaGLMpValues.mat
% LickMod=length(find(SuaGLMpValues(Condition1,8)<0.05)==1)/length(Condition1)



% cd(['T:\jan\HeadFixed Data\JK16\JK16_AATC2_lHPC_128_2017-11-13_11-26-21\spikes'])
% 
% Selection=[1,39];
% 
% [PsthSua,EventMatrix,EventType]=FENS_AATC_SuaPsth([1,2],1,4,50,1);
% 



%% Analysis sketches
clearvars -except TraceUpHPC TraceDownHPC          
cd('T:\jan\Collabo Data\PFCpaperPreProcessed')
load AATC_PFC_Sua_Psth_1ms

load('LickEvokedIndx.mat')
AATC_Sua_Psth=AATC_PFC_Sua_Psth;



%% Trace Down
Bin=25;
window=3000/Bin:4000/Bin;
baseline=1:1000/Bin;
time=-1+0.001*Bin:.001*Bin:4;

AATC_Sua_PsthBined=squeeze(mean(reshape(AATC_Sua_Psth,Bin,size(AATC_Sua_Psth,1)/Bin,size(AATC_Sua_Psth,2),size(AATC_Sua_Psth,3))));
%Subtract Baseline and Normalizse

bc_Psths=(AATC_Sua_PsthBined-nanmean(AATC_Sua_PsthBined(baseline,:,:)));%./nanmean(AATC_Sua_PsthBined(baseline,:,:))*100;;  %in herz

for i=1:size(bc_Psths,2)
   bc_Psths(:,i,:) = smoothdata(bc_Psths(:,i,:),'gaussian',25);
end

Tresh=-1;
BaseStd=squeeze(std(bc_Psths(baseline,:,:)));
EvR=squeeze(mean(bc_Psths(window,:,:)));
EvokedPeaks=squeeze(mean(bc_Psths(window,:,:)));
Evokedup(1:size(AATC_Sua_Psth,2))=0;
Evokedup(find(EvR(:,1)<BaseStd(:,1)*Tresh|EvR(:,2)<BaseStd(:,2)*Tresh))=1;

TraceDownPFC=Evokedup;

figure()
plotwindow=1:5000/Bin;
Condition1=find(LearnedCounter==1&Evokedup==1&LickDownPFC==0&LickUpPFC==0);
Condition=Condition1 

subplot(2,3,[1,2,4,5])

hold on
pp(1,:)=nanmean(bc_Psths(:,Condition,1)');
errBar=repmat(std(bc_Psths(:,Condition,2)')/sqrt(size(Condition,2)),2,1);
shadePlot2(time(plotwindow),pp(1,plotwindow),errBar(:,plotwindow),[0 0 1],.3)
hold on
pp(2,:)=nanmean(bc_Psths(:,Condition,2)');
errBar=repmat(std(bc_Psths(:,Condition,2)')/sqrt(size(Condition,2)),2,1);
shadePlot2(time(plotwindow),pp(2,plotwindow),errBar(:,plotwindow),[1 0 0],.3)
hold on
plot([0 0],[max(max(pp))+max(max(errBar)) min(min(pp))-min(min(errBar))],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
plot([2 2],[max(max(pp))+max(max(errBar)) min(min(pp))-min(min(errBar))],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
plot([3 3],[max(max(pp))+max(max(errBar)) min(min(pp))-min(min(errBar))],'LineWidth',2,'LineStyle',':','Color',[0 0 0])
box off
axis tight
xlabel('Time [s]')
ylabel('{\Delta} FiringRate [hz]')
set(gca,'FontSize',25);
set(gca,'LineWidth',5);


subplot(2,3,[3,6])
data{1}=EvokedPeaks(Condition,1);
data{2}=EvokedPeaks(Condition,2);
[p1,h]=ranksum(data{1},data{2})
y=[mean(data{1}),mean(data{2})]
e=[std(data{1})/sqrt(length(data{1})),...
std(data{2})/sqrt(length(data{2}))];
p=[ NaN p1;
    p1,NaN]
c=[0 0 1;1 0 0]
superbar(y,'E',e,'P',p,'BarFaceColor',c)
names={'CS+','CS-'}
set(gca,'Xtick',0:1:3)
set(gca,'xticklabel',names)
legend boxoff

% subplot(2,4,[8])
% NT(1)=length(find(Evokedup==1&neuronType==1&LearnedCounter==1))/length(find(neuronType==1&LearnedCounter==1))
% NT(2)=length(find(Evokedup==1&neuronType==2&LearnedCounter==1))/length(find(neuronType==2&LearnedCounter==1))
% bar(NT)
% 
% subplot(2,4,[8])
% NT(1)=length(find(Evokedup==1&neuronType==1&LearnedCounter==1))/length(find(Evokedup==1&LearnedCounter==1))
% NT(2)=length(find(Evokedup==1&neuronType==2&LearnedCounter==1))/length(find(Evokedup==1&LearnedCounter==1))
% NT(3)=length(find(Evokedup==1&(neuronType==0|neuronTypeGmm==3)&LearnedCounter==1))/length(find(Evokedup==1&LearnedCounter==1))
% pie(NT,[1 2 3])



clearvars data
Condition1=find(LearnedCounter==1&Evokedup==1&LickDownPFC==0&LickUpPFC==0);
Condition2=find(LearnedCounter==0&TrgDayCounter<3&Evokedup==1);

data{1}=EvokedPeaks(Condition1,1)';
data{2}=EvokedPeaks(Condition1,2)';
data{3}=EvokedPeaks(Condition2,1)';
data{4}=EvokedPeaks(Condition2,2)';

Reward=zeros(size([data{1},data{2},data{3},data{4}],2),1);
Reward([1:size(data{1},2),(size(data{1},2)*2)+1:size(data{1},2)*2+1+size(data{3},2)])=1;
Learned=zeros(size([data{1},data{2},data{3},data{4}],2),1);
Learned(1:size(data{1},2)+size(data{2},2))=1;
% [h,p]=anovan([data{1},data{2},data{3},data{4}],{Reward;Learned},'model','interaction')


AllPostCells=length(find(LearnedCounter==1))
PostTraceDown=length(find(LearnedCounter==1&Evokedup==1&LickDownPFC==0&LickUpPFC==0))
PerPostTracedown=PostTraceDown/AllPostCells
% figure()
% scatter(LickRateCounter(Evokedup==1),EvokedPeaks(Evokedup==1,1)-EvokedPeaks(Evokedup==1,2))
% lsline
% axis tight
% 
% [h,p]=corrcoef(LickRateCounter(Evokedup==1),EvokedPeaks(Evokedup==1,1)-EvokedPeaks(Evokedup==1,2))
% 
% load SuaGLMpValues.mat
% LickMod=length(find(SuaGLMpValues(Condition1,8)<0.05)==1)/length(Condition1)

%% Trace Up
Bin=25;
window=3000/Bin:4000/Bin;
baseline=1:1000/Bin;
time=-1+0.001*Bin:.001*Bin:4;

AATC_Sua_PsthBined=squeeze(mean(reshape(AATC_Sua_Psth,Bin,size(AATC_Sua_Psth,1)/Bin,size(AATC_Sua_Psth,2),size(AATC_Sua_Psth,3))));
%Subtract Baseline and Normalizse

bc_Psths=((AATC_Sua_PsthBined)-nanmean(AATC_Sua_PsthBined(baseline,:,:)));%./nanmean(AATC_Sua_PsthBined(baseline,:,:))*100;  %in herz

for i=1:size(bc_Psths,2)
   bc_Psths(:,i,:) = smoothdata(bc_Psths(:,i,:),'gaussian',25);
end

Tresh=1;
BaseStd=squeeze(std(bc_Psths(baseline,:,:)));
EvR=squeeze(mean((bc_Psths(window,:,:))));
EvokedPeaks=squeeze(mean(bc_Psths(window,:,:)));
Evokedup(1:size(AATC_Sua_Psth,2))=0;
Evokedup(find(EvR(:,1)>BaseStd(:,1)*Tresh|EvR(:,2)>BaseStd(:,2)*Tresh))=1;
MeanBaselineX=squeeze(mean(AATC_Sua_PsthBined(baseline,:,1)));

TraceUpPFC=Evokedup;

figure()

plotwindow=1:5000/Bin;
Condition1=find(LearnedCounter==1&Evokedup==1&LickDownPFC==0&LickUpPFC==0);
Condition=Condition1
mean(FRperiodsCounter(LearnedCounter==1&Evokedup==0,1))

subplot(2,3,[1,2,4,5])

hold on
pp(1,:)=nanmean(bc_Psths(:,Condition,1)');
errBar=repmat(std(bc_Psths(:,Condition,2)')/sqrt(size(Condition,2)),2,1);
shadePlot2(time(plotwindow),pp(1,plotwindow),errBar(:,plotwindow),[0 0 1],.3)
hold on
pp(2,:)=nanmean(bc_Psths(:,Condition,2)');
errBar=repmat(std(bc_Psths(:,Condition,2)')/sqrt(size(Condition,2)),2,1);
shadePlot2(time(plotwindow),pp(2,plotwindow),errBar(:,plotwindow),[1 0 0],.3)
hold on
plot([0 0],[max(max(pp))+max(max(errBar)) min(min(pp))-min(min(errBar))],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
plot([2 2],[max(max(pp))+max(max(errBar)) min(min(pp))-min(min(errBar))],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
plot([3 3],[max(max(pp))+max(max(errBar)) min(min(pp))-min(min(errBar))],'LineWidth',2,'LineStyle',':','Color',[0 0 0])
box off
axis tight
xlabel('Time [s]')
ylabel('{\Delta} FiringRate [hz]')
set(gca,'FontSize',25);
set(gca,'LineWidth',5);


subplot(2,3,[3,6])
data{1}=EvokedPeaks(Condition,1);
data{2}=EvokedPeaks(Condition,2);
[p1,h]=ranksum(data{1},data{2})
y=[mean(data{1}),mean(data{2})]
e=[std(data{1})/sqrt(length(data{1})),...
std(data{2})/sqrt(length(data{2}))];
p=[ NaN p1;
    p1,NaN]
c=[0 0 1;1 0 0]
superbar(y,'E',e,'P',p,'BarFaceColor',c)
names={'CS+','CS-'}
set(gca,'Xtick',0:1:3)
set(gca,'xticklabel',names)
legend boxoff


% subplot(2,4,[8])
% NT(1)=length(find(Evokedup==1&neuronType==1&LearnedCounter==1))/length(find(neuronType==1&LearnedCounter==1))
% NT(2)=length(find(Evokedup==1&neuronType==2&LearnedCounter==1))/length(find(neuronType==2&LearnedCounter==1))
% bar(NT)
% subplot(2,4,[8])
% NT(1)=length(find(Evokedup==1&neuronType==1&LearnedCounter==1))/length(find(Evokedup==1&LearnedCounter==1))
% NT(2)=length(find(Evokedup==1&neuronType==2&LearnedCounter==1))/length(find(Evokedup==1&LearnedCounter==1))
% NT(3)=length(find(Evokedup==1&(neuronType==0|neuronTypeGmm==3)&LearnedCounter==1))/length(find(Evokedup==1&LearnedCounter==1))
% pie(NT,[1 2 3])



clearvars data
Condition1=find(LearnedCounter==1&Evokedup==1&LickDownPFC==0&LickUpPFC==0);
Condition2=find(LearnedCounter==0&TrgDayCounter<3&Evokedup==1);

data{1}=EvokedPeaks(Condition1,1)';
data{2}=EvokedPeaks(Condition1,2)';
data{3}=EvokedPeaks(Condition2,1)';
data{4}=EvokedPeaks(Condition2,2)';

Reward=zeros(size([data{1},data{2},data{3},data{4}],2),1);
Reward([1:size(data{1},2),(size(data{1},2)*2)+1:size(data{1},2)*2+1+size(data{3},2)])=1;
Learned=zeros(size([data{1},data{2},data{3},data{4}],2),1);
Learned(1:size(data{1},2)+size(data{2},2))=1;
% [h,p]=anovan([data{1},data{2},data{3},data{4}],{Reward;Learned},'model','interaction')

AllPostCells=length(find(LearnedCounter==1))
PostTraceUp=length(find(LearnedCounter==1&Evokedup==1&LickDownPFC==0&LickUpPFC==0))
PerPostTraceUp=PostTraceUp/AllPostCells
% figure()
% scatter(LickRateCounter(Evokedup==1),EvokedPeaks(Evokedup==1,1)-EvokedPeaks(Evokedup==1,2))
% lsline
% axis tight
% 
% [h,p]=corrcoef(LickRateCounter(Evokedup==1),EvokedPeaks(Evokedup==1,1)-EvokedPeaks(Evokedup==1,2))
% 
% 
% %%
% 
% load SuaGLMpValues.mat
% LickMod=length(find(SuaGLMpValues(Condition1,8)<0.05)==1)/length(Condition1)


% 
% cd(['T:\jan\HeadFixed Data\JK16\JK16_AATC2_lHPC_128_2017-11-13_11-26-21\spikes'])
% 
% Selection=[1,39];
% 
% [PsthSua,EventMatrix,EventType]=FENS_AATC_SuaPsth([1,2],1,4,50,1);


