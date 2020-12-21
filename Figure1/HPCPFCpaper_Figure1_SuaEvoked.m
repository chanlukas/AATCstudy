

%% Analysis sketches
clear all      
cd('T:\jan\Collabo Data\HPCpaperPreProcessed')
load AATC_Sua_Psth_1ms


%% Pre Process Single Cell Psths
Bin=25;
window=1000/Bin:1350/Bin;
baseline=1:1000/Bin;

time=-1+0.001*Bin:.001*Bin:4;

AATC_Sua_PsthBined=squeeze(mean(reshape(AATC_Sua_Psth,Bin,size(AATC_Sua_Psth,1)/Bin,size(AATC_Sua_Psth,2),size(AATC_Sua_Psth,3))));
%Subtract Baseline and Normalizse

bc_Psths=(AATC_Sua_PsthBined-nanmean(AATC_Sua_PsthBined(baseline,:,:)));%./nanmean(AATC_Sua_PsthBined(baseline,:,:));  %in herz
for i=1:size(bc_Psths,2)
   bc_Psths(:,i,:) = smoothdata(bc_Psths(:,i,:),'gaussian',20);
end

BaseStd=squeeze(std(bc_Psths(baseline,:,:)));
EvR=squeeze(mean(bc_Psths(window,:,:)));
EvokedPeaks=squeeze(max(bc_Psths(window,:,:)));
Evokedup(1:size(AATC_Sua_Psth,2))=0;
Evokedup(find(EvR(:,1)>BaseStd(:,1)*3|EvR(:,2)>BaseStd(:,2)*3))=1;
length(find(Evokedup==1))/size(Evokedup,2)
y=find(Evokedup==1)

figure()
Condition=find(LearnedCounter==1&TrgDayCounter>3);


subplot(1,4,[1 2])
hold on
pp(1,:)=nanmean(bc_Psths(:,Condition,1)');
errBar=repmat(std(bc_Psths(:,Condition,1)')/sqrt(size(Condition,2)),2,1);
shadePlot2(time,pp(1,:),errBar(:,:),[0 0 1],.3)
hold on
pp(2,:)=nanmean(bc_Psths(:,Condition,2)');
errBar=repmat(std(bc_Psths(:,Condition,2)')/sqrt(size(Condition,2)),2,1);
shadePlot2(time,pp(2,:),errBar(:,:),[1 0 0],.3)
hold on
plot([0 0],[-.75 1.25],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
plot([2 2],[-.75 1.25],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
plot([3 3],[-.75 1.25],'LineWidth',2,'LineStyle',':','Color',[0 0 0])

% plot([3 3],[max(max(pp))+max(max(errBar)) min(min(pp))-min(min(errBar))],'LineWidth',2,'LineStyle',':','Color',[0 0 0])
box off
% axis ([-1 4 -.75 1.15])
xlabel('Time [s]')
ylabel('{\Delta} FiringRate [hz]')
set(gca,'FontSize',25);
set(gca,'LineWidth',5);


subplot(1,4,[3])
data{1}=EvR(Condition,1);
data{2}=EvR(Condition,2);
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


subplot(1,4,[4])
data{1}=squeeze(mean(bc_Psths(3000/Bin:4000/Bin,Condition,1)));
data{2}=squeeze(mean(bc_Psths(3000/Bin:4000/Bin,Condition,2)));
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
[p1,h]=signrank(data{1})
[p1,h]=signrank(data{2})


data{1}=EvokedPeaks((LearnedCounter==0&TrgDayCounter<3),1)';
data{2}=EvokedPeaks((LearnedCounter==0&TrgDayCounter<3),2)';
data{3}=EvokedPeaks((LearnedCounter==1&TrgDayCounter>3),1)';
data{4}=EvokedPeaks((LearnedCounter==1&TrgDayCounter>3),2)';
Reward=zeros(size([data{1},data{2},data{3},data{4}],2),1);
Reward([1:size(data{1},2),(size(data{1},2)*2)+1:size(data{1},2)*2+1+size(data{3},2)])=1;
Learned=zeros(size([data{1},data{2},data{3},data{4}],2),1);
Learned(1:size(data{1},2)+size(data{2},2))=1;
[h,p]=anovan([data{1},data{2},data{3},data{4}],{Reward;Learned},'model','interaction')

%% Analysis sketches
clear all           
cd('T:\jan\Collabo Data\PFCpaperPreProcessed')
load AATC_PFC_Sua_Psth_1ms


AATC_Sua_Psth=AATC_PFC_Sua_Psth;
%% Pre Process Single Cell Psths
Bin=25;
window=1000/Bin:1350/Bin;
baseline=1:1000/Bin;

time=-1+0.001*Bin:.001*Bin:4;

AATC_Sua_PsthBined=squeeze(mean(reshape(AATC_Sua_Psth,Bin,size(AATC_Sua_Psth,1)/Bin,size(AATC_Sua_Psth,2),size(AATC_Sua_Psth,3))));
%Subtract Baseline and Normalizse

bc_Psths=(AATC_Sua_PsthBined-nanmean(AATC_Sua_PsthBined(baseline,:,:)));%./nanmean(AATC_Sua_PsthBined(baseline,:,:));;  %in herz

for i=1:size(bc_Psths,2)
   bc_Psths(:,i,:) = smoothdata(bc_Psths(:,i,:),'gaussian',20);
end

BaseStd=squeeze(std(bc_Psths(baseline,:,:)));
EvR=squeeze(mean(bc_Psths(window,:,:)));
EvokedPeaks=squeeze(max(bc_Psths(window,:,:)));
Evokedup(1:size(AATC_Sua_Psth,2))=0;
Evokedup(find(EvR(:,1)>BaseStd(:,1)*3|EvR(:,2)>BaseStd(:,2)*3))=1;
length(find(Evokedup==1))/size(Evokedup,2)
y=find(Evokedup==1)

figure()
Condition=find(LearnedCounter==1&TrgDayCounter>3);


subplot(1,4,[1 2])
hold on
pp(1,:)=nanmean(bc_Psths(:,Condition,1)');
errBar=repmat(std(bc_Psths(:,Condition,1)')/sqrt(size(Condition,2)),2,1);
shadePlot2(time,pp(1,:),errBar(:,:),[0 0 1],.3)
hold on
pp(2,:)=nanmean(bc_Psths(:,Condition,2)');
errBar=repmat(std(bc_Psths(:,Condition,2)')/sqrt(size(Condition,2)),2,1);
shadePlot2(time,pp(2,:),errBar(:,:),[1 0 0],.3)
hold on
plot([0 0],[-.25 2.6],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
plot([2 2],[-.25 2.6],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
plot([3 3],[-.25 2.6],'LineWidth',2,'LineStyle',':','Color',[0 0 0])
% plot([3 3],[max(max(pp))+max(max(errBar)) min(min(pp))-min(min(errBar))],'LineWidth',2,'LineStyle',':','Color',[0 0 0])
box off
% axis ([-1 4 -.25 2.6])
xlabel('Time [s]')
ylabel('{\Delta} FiringRate [hz]')
set(gca,'FontSize',25);
set(gca,'LineWidth',5);



subplot(1,4,[3])
data{1}=EvR(Condition,1);
data{2}=EvR(Condition,2);
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


subplot(1,4,[4])

%Trace Comparison

data{1}=squeeze(mean(bc_Psths(3000/Bin:4000/Bin,Condition,1)));
data{2}=squeeze(mean(bc_Psths(3000/Bin:4000/Bin,Condition,2)));
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
[p1,h]=signrank(data{1})
[p1,h]=signrank(data{2})

%%
data{1}=EvokedPeaks((LearnedCounter==0&TrgDayCounter<3),1)';
data{2}=EvokedPeaks((LearnedCounter==0&TrgDayCounter<3),2)';
data{3}=EvokedPeaks((LearnedCounter==1&TrgDayCounter>3),1)';
data{4}=EvokedPeaks((LearnedCounter==1&TrgDayCounter>3),2)';
Reward=zeros(size([data{1},data{2},data{3},data{4}],2),1);
Reward([1:size(data{1},2),(size(data{1},2)*2)+1:size(data{1},2)*2+1+size(data{3},2)])=1;
Learned=zeros(size([data{1},data{2},data{3},data{4}],2),1);
Learned(1:size(data{1},2)+size(data{2},2))=1;
[h,p]=anovan([data{1},data{2},data{3},data{4}],{Reward;Learned},'model','interaction')
