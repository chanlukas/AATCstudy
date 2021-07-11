 

%% Analysis sketches
clear all      
cd('/Users/Jan/Documents/PFC Paper/AATCdata/PreprocessedDATA_PFCpaper')
load AATC_Sua_Psth_1ms


%% Pre Process Single Cell Psths
Bin=25;
window=1000/Bin:1350/Bin;
baseline=1:5000/Bin;



time=-5+0.001*Bin:.001*Bin:25;

AATC_Sua_PsthBined=squeeze(nanmean(reshape(AATC_Sua_Psth,Bin,size(AATC_Sua_Psth,1)/Bin,size(AATC_Sua_Psth,2),size(AATC_Sua_Psth,3))));
%Subtract Baseline and Normalizse

BaselineSTD=nanstd(AATC_Sua_PsthBined(baseline,:,:));
BaselineSTD(find(BaselineSTD==0))=NaN;
BaselineMean=nanmean(AATC_Sua_PsthBined(baseline,:,:));


bc_Psths=(AATC_Sua_PsthBined-BaselineMean)./BaselineSTD; %in herz
for i=1:size(bc_Psths,2)
   bc_Psths(:,i,:) = smoothdata(bc_Psths(:,i,:),'gaussian',20);
end



figure()
Condition=find(LearnedCounter==1&TrgDayCounter>3);


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
axis tight
xlabel('Time [s]')
ylabel('{\Delta} FiringRate [hz]')
set(gca,'FontSize',25);
set(gca,'LineWidth',5);



%% Analysis sketches
clear all           
load AATC_PFC_Sua_Psth_1ms


AATC_Sua_Psth=AATC_PFC_Sua_Psth;
%% Pre Process Single Cell Psths
Bin=25;
window=1000/Bin:1350/Bin;
baseline=1:1000/Bin;


time=-5+0.001*Bin:.001*Bin:25;

AATC_Sua_PsthBined=squeeze(nanmean(reshape(AATC_Sua_Psth,Bin,size(AATC_Sua_Psth,1)/Bin,size(AATC_Sua_Psth,2),size(AATC_Sua_Psth,3))));
%Subtract Baseline and Normalizse

BaselineSTD=nanstd(AATC_Sua_PsthBined(baseline,:,:));
BaselineSTD(find(BaselineSTD==0))=NaN;
BaselineMean=nanmean(AATC_Sua_PsthBined(baseline,:,:));


bc_Psths=(AATC_Sua_PsthBined-BaselineMean)./BaselineSTD; %in herz

for i=1:size(bc_Psths,2)
   bc_Psths(:,i,:) = smoothdata(bc_Psths(:,i,:),'gaussian',25);
end


figure()
Condition=find(LearnedCounter==1&TrgDayCounter>3);


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
axis tight
xlabel('Time [s]')
ylabel('{\Delta} FiringRate [hz]')
set(gca,'FontSize',25);
set(gca,'LineWidth',5);



%%
clear all      
cd('/Users/Jan/Documents/PFC Paper/AATCdata/PreprocessedDATA_PFCpaper')
load AATC_Sua_Psth_1ms
load('LickEvokedIndx.mat')

Bin=25;
AATC_Sua_Psth=AATC_Sua_Psth(1:end,:,:);

AATC_Sua_Psth=squeeze(mean(reshape(AATC_Sua_Psth,Bin,size(AATC_Sua_Psth,1)/Bin,size(AATC_Sua_Psth,2),size(AATC_Sua_Psth,3))));
bc_Psths=(AATC_Sua_Psth);%-mean(AATC_Sua_Psth(1:5000/Bin,:,:)))./std(AATC_Sua_Psth(1:5000/Bin,:,:));


for i=1:size(bc_Psths,2)
   bc_Psths(:,i,:) = smoothdata(bc_Psths(:,i,:),'gaussian',25);
end


counter=1;
c=1;
for a=unique(AnimalCounter)
for s=1:size(unique(SessionCounter(AnimalCounter==a)),2);
 
PsthTrace=bc_Psths(:,AnimalCounter==a&SessionCounter==s&LickUpHPC==0&LickDownHPC==0,:);

if size(find(TaskModCounter(SessionCounter==s&AnimalCounter==a)==1),2)>2 %include only sessions with more then 10 task responsive cells
cellnum=size(find(TaskModCounter(SessionCounter==s&AnimalCounter==a)==1),2);

    %for every Bin, calculate distance between population vector in
    %Rewarded vs Unrewarded Contions
for i=1:size(AATC_Sua_Psth,1)
    in=squeeze(PsthTrace(i,:,:))';
out = squareform(pdist(in));
Distance(i,counter)=out(1,2)./cellnum;
end


if max(LearnedCounter(SessionCounter==s&AnimalCounter==a))==0
Late(counter)=0;
else
Late(counter)=1;
if length(find(Late==1))==12
 a
 s
end
end
counter=counter+1;
incl(c)=1;
end
c=c+1;
end
end

%%
PreExample=bc_Psths(:,AnimalCounter==7&LearnedCounter==0&SessionCounter==1,:);
PostExample=bc_Psths(:,AnimalCounter==1&LearnedCounter==1&SessionCounter==5,:);

figure()

base=(1:(1000/Bin))
Tone=((1000/Bin):(3000/Bin));
Trace=((3000/Bin):(4000/Bin));
Reward=((4000/Bin):(5000/Bin));
time=-5+Bin/1000:Bin/1000:25;

figure()
hold on
    DistanceBC=(Distance-nanmean(Distance(1:(5000/Bin)-1,:)))./nanstd(Distance(1:(5000/Bin)-1,:));
     pp(1,:)=nanmean(DistanceBC(:,Late==1),2);
    errBar1=repmat((std(DistanceBC(:,Late==1)')/sqrt(size(find(Late==1),2))),2,1);
    shadePlot2(time,pp(1,:),errBar1,[0 0 0],.5,'-',3)
    hold on

plot([0 0],[min(min(pp))-min(min(errBar1)) max(max(pp))+max(max(errBar1))],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
plot([2 2],[min(min(pp))-min(min(errBar1)) max(max(pp))+max(max(errBar1))],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
plot([3 3],[min(min(pp))-min(min(errBar1)) max(max(pp))+max(max(errBar1))],'LineWidth',2,'LineStyle',':','Color',[0 0 0])
box off
axis tight
ylabel('z-Score distance')
xlabel('Time [s]')
set(gca,'LineWidth',5)
set(gca,'FontSize',25)

[p,h]=signrank(nanmean(DistanceBC(1120:1200,Late==1)))

figure()
time=-1+Bin/1000:Bin/1000:4;
clearvars DistanceBC pp errBar1

hold on
    DistanceBC=(Distance-nanmean(Distance(1:(5000/Bin)-1,:)))./nanstd(Distance(1:(5000/Bin)-1,:));
     pp(1,:)=nanmean(DistanceBC(161:360,Late==1),2);
    errBar1=repmat((std(DistanceBC(161:360,Late==1)')/sqrt(size(find(Late==1),2))),2,1);
    shadePlot2(time,pp(1,:),errBar1,[0 0 0],.5,'-',3)
    hold on

plot([0 0],[min(min(pp))-min(min(errBar1)) max(max(pp))+max(max(errBar1))],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
plot([2 2],[min(min(pp))-min(min(errBar1)) max(max(pp))+max(max(errBar1))],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
plot([3 3],[min(min(pp))-min(min(errBar1)) max(max(pp))+max(max(errBar1))],'LineWidth',2,'LineStyle',':','Color',[0 0 0])
box off
axis tight
ylabel('z-Score distance')
xlabel('Time [s]')
set(gca,'LineWidth',5)
set(gca,'FontSize',25)


%%
clear all           
load('LickEvokedIndx.mat')
load AATC_PFC_Sua_Psth_1ms


AATC_Sua_Psth=AATC_PFC_Sua_Psth;
AATC_Sua_Psth=AATC_Sua_Psth(1:end,:,:);

Bin=25;
AATC_Sua_Psth=squeeze(mean(reshape(AATC_Sua_Psth,Bin,size(AATC_Sua_Psth,1)/Bin,size(AATC_Sua_Psth,2),size(AATC_Sua_Psth,3))));

% bc_Psths=AATC_Sua_Psth-mean(AATC_Sua_Psth(1:1000/Bin,:,:));

bc_Psths=(AATC_Sua_Psth);%-mean(AATC_Sua_Psth(1:5000/Bin,:,:)))./std(AATC_Sua_Psth(1:5000/Bin,:,:));


for i=1:size(bc_Psths,2)
   bc_Psths(:,i,:) = smoothdata(bc_Psths(:,i,:),'gaussian',25);
end


counter=1;
c=1;
for a=unique(AnimalCounter)
for s=1:size(unique(SessionCounter(AnimalCounter==a)),2);
 
PsthTrace=bc_Psths(:,AnimalCounter==a&SessionCounter==s&LickUpPFC==0&LickDownPFC==0,:);

if size(find(TaskModCounter(SessionCounter==s&AnimalCounter==a)==1),2)>2 %include only sessions with more then 10 task responsive cells
cellnum=size(find(TaskModCounter(SessionCounter==s&AnimalCounter==a)==1),2);

    %for every Bin, calculate distance between population vector in
    %Rewarded vs Unrewarded Contions
for i=1:size(AATC_Sua_Psth,1)
    in=squeeze(PsthTrace(i,:,:))';
out = squareform(pdist(in));
Distance(i,counter)=out(1,2)./cellnum;
end


if max(LearnedCounter(SessionCounter==s&AnimalCounter==a))==0
Late(counter)=0;
else
Late(counter)=1;
if length(find(Late==1))==12
 a
 s
end
end
counter=counter+1;
incl(c)=1;
end
c=c+1;
end
end

%%
PostExample=bc_Psths(:,find(AnimalCounter==16&LearnedCounter==1&SessionCounter==4),:);

figure()

base=(1:(1000/Bin))
Tone=((1000/Bin):(3000/Bin));
Trace=((3000/Bin):(4000/Bin));
Reward=((4000/Bin):(5000/Bin));
time=-5+Bin/1000:Bin/1000:25;

figure()
hold on
    DistanceBC=(Distance-nanmean(Distance(1:(5000/Bin)-1,:)))./nanstd(Distance(1:(5000/Bin)-1,:));
     pp(1,:)=nanmean(DistanceBC(:,Late==1),2);
    errBar1=repmat((std(DistanceBC(:,Late==1)')/sqrt(size(find(Late==1),2))),2,1);
    shadePlot2(time,pp(1,:),errBar1,[0 0 0],.5,'-',3)
    hold on

plot([0 0],[min(min(pp))-min(min(errBar1)) max(max(pp))+max(max(errBar1))],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
plot([2 2],[min(min(pp))-min(min(errBar1)) max(max(pp))+max(max(errBar1))],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
plot([3 3],[min(min(pp))-min(min(errBar1)) max(max(pp))+max(max(errBar1))],'LineWidth',2,'LineStyle',':','Color',[0 0 0])
box off
axis tight
ylabel('z-Score distance')
xlabel('Time [s]')
set(gca,'LineWidth',5)
set(gca,'FontSize',25)

[p,h]=signrank(nanmean(DistanceBC(1120:1200,Late==1)))


figure()
time=-1+Bin/1000:Bin/1000:4;
clearvars DistanceBC pp errBar1
hold on
    DistanceBC=(Distance-nanmean(Distance(1:(5000/Bin)-1,:)))./nanstd(Distance(1:(5000/Bin)-1,:));
     pp(1,:)=nanmean(DistanceBC(161:360,Late==1),2);
    errBar1=repmat((std(DistanceBC(161:360,Late==1)')/sqrt(size(find(Late==1),2))),2,1);
    shadePlot2(time,pp(1,:),errBar1,[0 0 0],.5,'-',3)
    hold on

plot([0 0],[min(min(pp))-min(min(errBar1)) max(max(pp))+max(max(errBar1))],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
plot([2 2],[min(min(pp))-min(min(errBar1)) max(max(pp))+max(max(errBar1))],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
plot([3 3],[min(min(pp))-min(min(errBar1)) max(max(pp))+max(max(errBar1))],'LineWidth',2,'LineStyle',':','Color',[0 0 0])
box off
axis tight
ylabel('z-Score distance')
xlabel('Time [s]')
set(gca,'LineWidth',5)
set(gca,'FontSize',25)