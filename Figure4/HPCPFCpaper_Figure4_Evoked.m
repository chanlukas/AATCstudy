

%% Analysis sketches
clear all   

cd('T:\jan\HeadFixed Data\AATCR2_HPC128-PFC128\AATCR2_Processed_Data')
load AATCR2_Sua_Psth_1ms
AATCR2_Sua_Psth(:,find(ismember(AnimalCounter,[2,4])),[1,2,3])=AATCR2_Sua_Psth(:,find(ismember(AnimalCounter,[2,4])),[3,2,1]);

AATC_Sua_Psth=AATCR2_Sua_Psth(1:5000,:,:);


%% Pre Process Single Cell Psths
Bin=25;
window=1000/Bin:1300/Bin;
baseline=1:1000/Bin;

time=-1+0.001*Bin:.001*Bin:4;

AATC_Sua_PsthBined=squeeze(mean(reshape(AATC_Sua_Psth,Bin,size(AATC_Sua_Psth,1)/Bin,size(AATC_Sua_Psth,2),size(AATC_Sua_Psth,3))));

bc_Psths=(AATC_Sua_PsthBined-nanmean(AATC_Sua_PsthBined(baseline,:,:)));  %in herz

for i=1:size(bc_Psths,2)
   bc_Psths(:,i,:) = smoothdata(bc_Psths(:,i,:),'gaussian',20);
end

BaseStd=squeeze(std(bc_Psths(baseline,:,:)));
EvR=squeeze(max(bc_Psths(window,:,:)));
EvokedPeaks=squeeze(mean(bc_Psths(window,:,:)));
Evokedup(1:size(AATC_Sua_Psth,2))=0;
Evokedup(find(EvR(:,1)>BaseStd(:,1)*2|EvR(:,2)>BaseStd(:,2)*2|EvR(:,3)>BaseStd(:,3)*2))=1;


plotwindow=1:5000/Bin;
Condition1=find(Learned2Counter==1&LearnedCounter==1&ShankS==1);
Condition=Condition1 

figure()

subplot(2,4,[1,2,5,6])
hold on
pp(1,:)=nanmean(bc_Psths(:,Condition,1)');
errBar=repmat(std(bc_Psths(:,Condition,2)')/sqrt(size(Condition,2)),2,1);
shadePlot2(time(plotwindow),pp(1,plotwindow),errBar(:,plotwindow),[0 0 1],.3)
hold on
pp(2,:)=nanmean(bc_Psths(:,Condition,2)');
errBar=repmat(std(bc_Psths(:,Condition,2)')/sqrt(size(Condition,2)),2,1);
shadePlot2(time(plotwindow),pp(2,plotwindow),errBar(:,plotwindow),[1 0 0],.3)
hold on
pp(3,:)=nanmean(bc_Psths(:,Condition,3)');
errBar=repmat(std(bc_Psths(:,Condition,3)')/sqrt(size(Condition,2)),2,1);
shadePlot2(time(plotwindow),pp(3,plotwindow),errBar(:,plotwindow),[0 0 .7],.3)
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

subplot(2,4,[3,7])
data{1}=EvokedPeaks(Condition,1);
data{2}=EvokedPeaks(Condition,2);
data{3}=EvokedPeaks(Condition,3);
[p1,h]=ranksum(data{1},data{2});
[p2,h]=ranksum(data{2},data{3});
[p3,h]=ranksum(data{1},data{3})

y=[mean(data{1}),mean(data{2}),mean(data{3})]

e=[std(data{1})/sqrt(length(data{1})),...
 std(data{2})/sqrt(length(data{2})),...
std(data{3})/sqrt(length(data{3}))]

p=[ NaN p1 NaN;...
    p1 NaN p2;...
    NaN p2 NaN]

c=[1 0 0;0 0 1;1 0 0]
superbar(y,'E',e,'P',p,'BarFaceColor',c)
box off


subplot(2,4,[4])

RewardSoundSpecific=find(ShankS'==1&Learned2Counter'==1&((EvR(:,1)>BaseStd(:,1)*2&EvR(:,3)>BaseStd(:,3)*2)&EvR(:,2)<BaseStd(:,2)*2));

SoundSpecific1=find(ShankS'==1&Learned2Counter'==1&((EvR(:,1)>BaseStd(:,1)*2&EvR(:,3)<BaseStd(:,3)*2&EvR(:,2)<BaseStd(:,2)*2)));
SoundSpecific2=find(ShankS'==1&Learned2Counter'==1&((EvR(:,2)>BaseStd(:,2)*2&EvR(:,1)<BaseStd(:,1)*2&EvR(:,3)<BaseStd(:,3)*2)));
SoundSpecific3=find(ShankS'==1&Learned2Counter'==1&((EvR(:,3)>BaseStd(:,3)*2&EvR(:,1)<BaseStd(:,1)*2&EvR(:,2)<BaseStd(:,2)*2)));

SoundUnSpecific=find(ShankS'==1&Learned2Counter'==1&(EvR(:,1)>BaseStd(:,1)*2&EvR(:,2)>BaseStd(:,2)*2&EvR(:,3)>BaseStd(:,3)*2));

RewardUnRewardCombi1=find(ShankS'==1&Learned2Counter'==1&((EvR(:,1)>BaseStd(:,1)*2&EvR(:,3)<BaseStd(:,3)*2&EvR(:,2)>BaseStd(:,2)*2)));
RewardUnRewardCombi2=find(ShankS'==1&Learned2Counter'==1&((EvR(:,1)<BaseStd(:,1)*2&EvR(:,3)>BaseStd(:,3)*2&EvR(:,2)>BaseStd(:,2)*2)));

RT(1)=length([SoundSpecific1;SoundSpecific2;SoundSpecific3])
RT(2)=length(RewardSoundSpecific)
RT(3)=length([RewardUnRewardCombi1;RewardUnRewardCombi2;SoundUnSpecific])
pie(RT,[1 2 3])

subplot(2,4,[8])
NT(1)=length(find(Evokedup==1&neuronType==1&Learned2Counter==1&LearnedCounter==1&ShankS==1))/length(find(Evokedup==1&Learned2Counter==1&LearnedCounter==1&ShankS==1))
NT(2)=length(find(Evokedup==1&neuronType==2&Learned2Counter==1&LearnedCounter==1&ShankS==1))/length(find(Evokedup==1&Learned2Counter==1&LearnedCounter==1&ShankS==1))
NT(3)=length(find(Evokedup==1&neuronType==0&Learned2Counter==1&LearnedCounter==1&ShankS==1))/length(find(Evokedup==1&Learned2Counter==1&LearnedCounter==1&ShankS==1))
pie(NT,[1 2 3])


clearvars RT NT RewardSoundSpecific SoundSpecific1 SoundSpecific2 SoundSpecific3 SoundUnSpecific RewardUnRewardCombi1 RewardUnRewardCombi2
figure()

plotwindow=1:5000/Bin;
Condition1=find(Learned2Counter==1&LearnedCounter==1&ShankS==2);
Condition=Condition1 

subplot(2,4,[1,2,5,6])
hold on
pp(1,:)=nanmean(bc_Psths(:,Condition,1)');
errBar=repmat(std(bc_Psths(:,Condition,2)')/sqrt(size(Condition,2)),2,1);
shadePlot2(time(plotwindow),pp(1,plotwindow),errBar(:,plotwindow),[0 0 1],.3)
hold on
pp(2,:)=nanmean(bc_Psths(:,Condition,2)');
errBar=repmat(std(bc_Psths(:,Condition,2)')/sqrt(size(Condition,2)),2,1);
shadePlot2(time(plotwindow),pp(2,plotwindow),errBar(:,plotwindow),[1 0 0],.3)
hold on
pp(3,:)=nanmean(bc_Psths(:,Condition,3)');
errBar=repmat(std(bc_Psths(:,Condition,3)')/sqrt(size(Condition,2)),2,1);
shadePlot2(time(plotwindow),pp(3,plotwindow),errBar(:,plotwindow),[0 0 .7],.3)
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

subplot(2,4,[3,7])
data{1}=EvokedPeaks(Condition,1);
data{2}=EvokedPeaks(Condition,2);
data{3}=EvokedPeaks(Condition,3);
[p1,h]=ranksum(data{1},data{2});
[p2,h]=ranksum(data{2},data{3});
[p3,h]=ranksum(data{1},data{3})


y=[mean(data{1}),mean(data{2}),mean(data{3})]

e=[std(data{1})/sqrt(length(data{1})),...
 std(data{2})/sqrt(length(data{2})),...
std(data{3})/sqrt(length(data{3}))]

p=[ NaN p1 NaN;...
    p1 NaN p2;...
    NaN p2 NaN]

c=[1 0 0;0 0 1;1 0 0]
superbar(y,'E',e,'P',p,'BarFaceColor',c)
box off

subplot(2,4,[4])

RewardSoundSpecific=find(ShankS'==2&Learned2Counter'==1&((EvR(:,1)>BaseStd(:,1)*2&EvR(:,3)>BaseStd(:,3)*2)&EvR(:,2)<BaseStd(:,2)*2));

SoundSpecific1=find(ShankS'==2&Learned2Counter'==1&((EvR(:,1)>BaseStd(:,1)*2&EvR(:,3)<BaseStd(:,3)*2&EvR(:,2)<BaseStd(:,2)*2)));
SoundSpecific2=find(ShankS'==2&Learned2Counter'==1&((EvR(:,2)>BaseStd(:,2)*2&EvR(:,1)<BaseStd(:,1)*2&EvR(:,3)<BaseStd(:,3)*2)));
SoundSpecific3=find(ShankS'==2&Learned2Counter'==1&((EvR(:,3)>BaseStd(:,3)*2&EvR(:,1)<BaseStd(:,1)*2&EvR(:,2)<BaseStd(:,2)*2)));

SoundUnSpecific=find(ShankS'==2&Learned2Counter'==1&(EvR(:,1)>BaseStd(:,1)*2&EvR(:,2)>BaseStd(:,2)*2&EvR(:,3)>BaseStd(:,3)*2));

RewardUnRewardCombi1=find(ShankS'==2&Learned2Counter'==1&((EvR(:,1)>BaseStd(:,1)*2&EvR(:,3)<BaseStd(:,3)*2&EvR(:,2)>BaseStd(:,2)*2)));
RewardUnRewardCombi2=find(ShankS'==2&Learned2Counter'==1&((EvR(:,1)<BaseStd(:,1)*2&EvR(:,3)>BaseStd(:,3)*2&EvR(:,2)>BaseStd(:,2)*2)));

RT(1)=length([SoundSpecific1;SoundSpecific2;SoundSpecific3])
RT(2)=length(RewardSoundSpecific)
RT(3)=length([RewardUnRewardCombi1;RewardUnRewardCombi2;SoundUnSpecific])

pie(RT,[1 2 3])

subplot(2,4,[8])
NT(1)=length(find(Evokedup==1&neuronType==1&Learned2Counter==1&LearnedCounter==1&ShankS==2))/length(find(Evokedup==1&Learned2Counter==1&LearnedCounter==1&ShankS==2))
NT(2)=length(find(Evokedup==1&neuronType==2&Learned2Counter==1&LearnedCounter==1&ShankS==2))/length(find(Evokedup==1&Learned2Counter==1&LearnedCounter==1&ShankS==2))
NT(3)=length(find(Evokedup==1&neuronType==0&ShankS==2&Learned2Counter==1&LearnedCounter==1))/length(find(Evokedup==1&Learned2Counter==1&LearnedCounter==1&ShankS==2))
pie(NT,[1 2 3])


