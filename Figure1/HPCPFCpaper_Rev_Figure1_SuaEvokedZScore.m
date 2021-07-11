

%% Analysis sketches
clear all      
cd('/Users/Jan/Documents/PFC Paper/AATCdata/PreprocessedDATA_PFCpaper')
load AATC_Sua_Psth_1ms

[x]=max(FRperiodsCounter,[],2);
[y]=min(FRperiodsCounter,[],2);
z=x./y;
TF = isoutlier(z,'quartiles');
%% Pre Process Single Cell Psths
Bin=25;
window=1000/Bin:1350/Bin;
baseline=1:1000/Bin;

AATC_Sua_Psth=AATC_Sua_Psth(4001:10000,:,:);


time=-1+0.001*Bin:.001*Bin:5;

AATC_Sua_PsthBined=squeeze(nanmean(reshape(AATC_Sua_Psth,Bin,size(AATC_Sua_Psth,1)/Bin,size(AATC_Sua_Psth,2),size(AATC_Sua_Psth,3))));
%Subtract Baseline and Normalizse

BaselineSTD=nanstd(AATC_Sua_PsthBined(baseline,:,:));
BaselineSTD(find(BaselineSTD==0))=NaN;
BaselineMean=nanmean(AATC_Sua_PsthBined(baseline,:,:));


bc_Psths=(AATC_Sua_PsthBined-BaselineMean)./BaselineSTD; %in herz
for i=1:size(bc_Psths,2)
   bc_Psths(:,i,:) = smoothdata(bc_Psths(:,i,:),'gaussian',20);
end

BaseStd=squeeze(nanstd(bc_Psths(baseline,:,:)));
EvR=squeeze(nanmean(bc_Psths(window,:,:)));
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
errBar=repmat(nanstd(bc_Psths(:,Condition,1)')/sqrt(size(Condition,2)),2,1);
shadePlot2(time,pp(1,:),errBar(:,:),[0 0 1],.3)
hold on
pp(2,:)=nanmean(bc_Psths(:,Condition,2)');
errBar=repmat(nanstd(bc_Psths(:,Condition,2)')/sqrt(size(Condition,2)),2,1);
shadePlot2(time,pp(2,:),errBar(:,:),[1 0 0],.3)
hold on
plot([0 0],[min(min(pp))-min(min(errBar)) max(max(pp))+max(max(errBar))],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
plot([2 2],[min(min(pp))-min(min(errBar)) max(max(pp))+max(max(errBar))],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
plot([3 3],[min(min(pp))-min(min(errBar)) max(max(pp))+max(max(errBar))],'LineWidth',2,'LineStyle',':','Color',[0 0 0])

% plot([3 3],[max(max(pp))+max(max(errBar)) min(min(pp))-min(min(errBar))],'LineWidth',2,'LineStyle',':','Color',[0 0 0])
box off
axis ([-1 4 -.25 .3])
xlabel('Time [s]')
ylabel('{\Delta} FiringRate [hz]')
set(gca,'FontSize',25);
set(gca,'LineWidth',5);


subplot(1,4,[3])
data{1}=EvR(Condition,1);
data{2}=EvR(Condition,2);
[p1,h]=ranksum(data{1},data{2})
y=[nanmean(data{1}),nanmean(data{2})]
e=[nanstd(data{1})/sqrt(length(data{1})),...
nanstd(data{2})/sqrt(length(data{2}))];
p=[ NaN p1;
    p1,NaN]
c=[0 0 1;1 0 0]
superbar(y,'E',e,'P',p,'BarFaceColor',c)
names={'CS+','CS-'}
set(gca,'Xtick',0:1:3)
set(gca,'xticklabel',names)
% legend boxoff


subplot(1,4,[4])
data{1}=squeeze(nanmean(bc_Psths(3000/Bin:4000/Bin,Condition,1)));
data{2}=squeeze(nanmean(bc_Psths(3000/Bin:4000/Bin,Condition,2)));
[p1,h]=ranksum(data{1},data{2})
y=[nanmean(data{1}),nanmean(data{2})]
e=[nanstd(data{1})/sqrt(length(data{1})),...
nanstd(data{2})/sqrt(length(data{2}))];
p=[ NaN p1;
    p1,NaN]
c=[0 0 1;1 0 0]
superbar(y,'E',e,'P',p,'BarFaceColor',c)
names={'CS+','CS-'}
set(gca,'Xtick',0:1:3)
set(gca,'xticklabel',names)
% legend boxoff
[p1,h]=signrank(data{1})
[p1,h]=signrank(data{2})

%% Test All Sound responses for normality
[h,p] = kstest( [EvokedPeaks(:,1);EvokedPeaks(:,2)] ) 


%kruskalwallis to check for group differences
data1=EvokedPeaks((LearnedCounter==0&TrgDayCounter<3),1)';
data2=EvokedPeaks((LearnedCounter==0&TrgDayCounter<3),2)';
data3=EvokedPeaks((LearnedCounter==1&TrgDayCounter>3),1)';
data4=EvokedPeaks((LearnedCounter==1&TrgDayCounter>3),2)';
data=cat(2,data1,data2,data3,data4);
label1(1:length(data1))=1;
label2(1:length(data2))=2;
label3(1:length(data3))=3;
label4(1:length(data4))=4;
label=cat(2,label1,label2,label3,label4);
p = kruskalwallis(data,label)

%%Pre vs Post CS+ and CS-
[p1,h]=ranksum(data1,data3)
[p1,h]=ranksum(data2,data4)
[p1,h]=ranksum(data3,data4)

sessionsCounter(LearnedCounter==0&TrgDayCounter<3)
%anova
% data{1}=EvokedPeaks((LearnedCounter==0&TrgDayCounter<3),1)';
% data{2}=EvokedPeaks((LearnedCounter==0&TrgDayCounter<3),2)';
% data{3}=EvokedPeaks((LearnedCounter==1&TrgDayCounter>3),1)';
% data{4}=EvokedPeaks((LearnedCounter==1&TrgDayCounter>3),2)';
% Reward=zeros(size([data{1},data{2},data{3},data{4}],2),1);
% Reward([1:size(data{1},2),(size(data{1},2)*2)+1:size(data{1},2)*2+1+size(data{3},2)])=1;
% Learned=zeros(size([data{1},data{2},data{3},data{4}],2),1);
% Learned(1:size(data{1},2)+size(data{2},2))=1;
% [h,p]=anovan([data{1},data{2},data{3},data{4}],{Reward;Learned},'model','interaction')


% %mackskill
% 
% data1=EvokedPeaks((LearnedCounter==0&TrgDayCounter<3),1)';
% data2=EvokedPeaks((LearnedCounter==0&TrgDayCounter<3),2)';
% data3=EvokedPeaks((LearnedCounter==1&TrgDayCounter>3),1)';
% data4=EvokedPeaks((LearnedCounter==1&TrgDayCounter>3),2)';
% 
% data12=cat(1,data1,data2);
% data12b=NaN(2,length(data34))
% data12b(:,1:length(data12))=data12;
% data34=cat(1,data3,data4);
% data=cat(2,data12b,data34);
% [p stats]=mackskill(data',length(data34))
% 
% 
% 
% %Friedman CS+ vs CS-
% 
% data1=EvokedPeaks((LearnedCounter==0&TrgDayCounter<3),1)';
% data2=EvokedPeaks((LearnedCounter==0&TrgDayCounter<3),2)';
% data3=EvokedPeaks((LearnedCounter==1&TrgDayCounter>3),1)';
% data4=EvokedPeaks((LearnedCounter==1&TrgDayCounter>3),2)';
% 
% data12=cat(1,data1,data2);
% data34=cat(1,data3,data4);
% subsample=randsample(length(data34),length(data12))
% data34=data34(:,subsample);
% data=cat(2,data12,data34);
% [p stats]=friedman(data',length(data12))
% 
% %Friedman Pre vs Post
% 
% data1=EvokedPeaks((LearnedCounter==0&TrgDayCounter<3),1)';
% data2=EvokedPeaks((LearnedCounter==0&TrgDayCounter<3),2)';
% data3=EvokedPeaks((LearnedCounter==1&TrgDayCounter>3),1)';
% data4=EvokedPeaks((LearnedCounter==1&TrgDayCounter>3),2)';
% 
% data12=cat(2,data1,data2);
% data34=cat(2,data3,data4);
% subsample=randsample(length(data34),length(data12))
% data34=data34(:,subsample);
% data=cat(1,data12,data34);
% [p stats]=friedman(data,length(data12))

%%
%Reward Comparison

Condition2=find(TrgDayCounter<3&LearnedCounter==0&TF'==0);

figure()
subplot(1,4,[4])

data{1}=squeeze(nanmean(bc_Psths(4000/Bin:4500/Bin,Condition2,1)));
data{2}=squeeze(nanmean(bc_Psths(4000/Bin:4500/Bin,Condition,1)));
[p1,h]=ranksum(data{1},data{2})
y=[nanmean(data{1}),nanmean(data{2})]
e=[nanstd(data{1})/sqrt(length(data{1})),...
nanstd(data{2})/sqrt(length(data{2}))];
p=[ NaN p1;
    p1,NaN]
c=[0 0 1;1 0 0]
superbar(y,'E',e,'P',p,'BarFaceColor',c)
names={'CS+','CS-'}
set(gca,'Xtick',0:1:3)
set(gca,'xticklabel',names)
% legend boxoff
[p1,h]=signrank(data{1})
[p1,h]=signrank(data{2})

%% Analysis sketches
clear all           
load AATC_PFC_Sua_Psth_1ms

[x]=max(FRperiodsCounter,[],2);
[y]=min(FRperiodsCounter,[],2);
z=x./y;
TF = isoutlier(z,'quartiles');

AATC_Sua_Psth=AATC_PFC_Sua_Psth;
%% Pre Process Single Cell Psths
Bin=25;
window=1000/Bin:1350/Bin;
baseline=1:1000/Bin;

AATC_Sua_Psth=AATC_Sua_Psth(4001:10000,:,:);

time=-1+0.001*Bin:.001*Bin:5;

AATC_Sua_PsthBined=squeeze(nanmean(reshape(AATC_Sua_Psth,Bin,size(AATC_Sua_Psth,1)/Bin,size(AATC_Sua_Psth,2),size(AATC_Sua_Psth,3))));
%Subtract Baseline and Normalizse

BaselineSTD=nanstd(AATC_Sua_PsthBined(baseline,:,:));
BaselineSTD(find(BaselineSTD==0))=NaN;
BaselineMean=nanmean(AATC_Sua_PsthBined(baseline,:,:));


bc_Psths=(AATC_Sua_PsthBined-BaselineMean)./BaselineSTD; %in herz

for i=1:size(bc_Psths,2)
   bc_Psths(:,i,:) = smoothdata(bc_Psths(:,i,:),'gaussian',20);
end

BaseStd=squeeze(nanstd(bc_Psths(baseline,:,:)));
EvR=squeeze(nanmean(bc_Psths(window,:,:)));
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
errBar=repmat(nanstd(bc_Psths(:,Condition,1)')/sqrt(size(Condition,2)),2,1);
shadePlot2(time,pp(1,:),errBar(:,:),[0 0 1],.3)
hold on
pp(2,:)=nanmean(bc_Psths(:,Condition,2)');
errBar=repmat(nanstd(bc_Psths(:,Condition,2)')/sqrt(size(Condition,2)),2,1);
shadePlot2(time,pp(2,:),errBar(:,:),[1 0 0],.3)
hold on
plot([0 0],[min(min(pp))-min(min(errBar)) max(max(pp))+max(max(errBar))],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
plot([2 2],[min(min(pp))-min(min(errBar)) max(max(pp))+max(max(errBar))],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
plot([3 3],[min(min(pp))-min(min(errBar)) max(max(pp))+max(max(errBar))],'LineWidth',2,'LineStyle',':','Color',[0 0 0])
% plot([3 3],[max(max(pp))+max(max(errBar)) min(min(pp))-min(min(errBar))],'LineWidth',2,'LineStyle',':','Color',[0 0 0])
box off
axis ([-1 4 -.25 1.5])
xlabel('Time [s]')
ylabel('{\Delta} FiringRate [hz]')
set(gca,'FontSize',25);
set(gca,'LineWidth',5);



subplot(1,4,[3])
data{1}=EvR(Condition,1);
data{2}=EvR(Condition,2);
[p1,h]=ranksum(data{1},data{2})
y=[nanmean(data{1}),nanmean(data{2})]
e=[nanstd(data{1})/sqrt(length(data{1})),...
nanstd(data{2})/sqrt(length(data{2}))];
p=[ NaN p1;
    p1,NaN]
c=[0 0 1;1 0 0]
superbar(y,'E',e,'P',p,'BarFaceColor',c)
names={'CS+','CS-'}
set(gca,'Xtick',0:1:3)
set(gca,'xticklabel',names)
% legend boxoff


subplot(1,4,[4])

%Trace Comparison

data{1}=squeeze(nanmean(bc_Psths(3000/Bin:4000/Bin,Condition,1)));
data{2}=squeeze(nanmean(bc_Psths(3000/Bin:4000/Bin,Condition,2)));
[p1,h]=ranksum(data{1},data{2})
y=[nanmean(data{1}),nanmean(data{2})]
e=[nanstd(data{1})/sqrt(length(data{1})),...
nanstd(data{2})/sqrt(length(data{2}))];
p=[ NaN p1;
    p1,NaN]
c=[0 0 1;1 0 0]
superbar(y,'E',e,'P',p,'BarFaceColor',c)
names={'CS+','CS-'}
set(gca,'Xtick',0:1:3)
set(gca,'xticklabel',names)
% legend boxoff
[p1,h]=signrank(data{1})
[p1,h]=signrank(data{2})

%% Test All Sound responses for normality
[h,p] = kstest( [EvokedPeaks(:,1);EvokedPeaks(:,2)] ) 


%kruskalwallis to check for group differences
data1=EvokedPeaks((LearnedCounter==0&TrgDayCounter<3),1)';
data2=EvokedPeaks((LearnedCounter==0&TrgDayCounter<3),2)';
data3=EvokedPeaks((LearnedCounter==1&TrgDayCounter>3),1)';
data4=EvokedPeaks((LearnedCounter==1&TrgDayCounter>3),2)';
data=cat(2,data1,data2,data3,data4);
label1(1:length(data1))=1;
label2(1:length(data2))=2;
label3(1:length(data3))=3;
label4(1:length(data4))=4;
label=cat(2,label1,label2,label3,label4);
p = kruskalwallis(data,label)

%%Pre vs Post CS+ and CS-
[p1,h]=ranksum(data1,data3)
[p1,h]=ranksum(data2,data4)

%%
%Reward Comparison

Condition2=find(TrgDayCounter<3&LearnedCounter==0);

figure()
subplot(1,4,[4])

data{1}=squeeze(nanmean(bc_Psths(4000/Bin:4500/Bin,Condition2,1)));
data{2}=squeeze(nanmean(bc_Psths(4000/Bin:4500/Bin,Condition,1)));
[p1,h]=ranksum(data{1},data{2})
y=[nanmean(data{1}),nanmean(data{2})]
e=[nanstd(data{1})/sqrt(length(data{1})),...
nanstd(data{2})/sqrt(length(data{2}))];
p=[ NaN p1;
    p1,NaN]
c=[0 0 1;1 0 0]
superbar(y,'E',e,'P',p,'BarFaceColor',c)
names={'CS+','CS-'}
set(gca,'Xtick',0:1:3)
set(gca,'xticklabel',names)
% legend boxoff
[p1,h]=signrank(data{1})
[p1,h]=signrank(data{2})
