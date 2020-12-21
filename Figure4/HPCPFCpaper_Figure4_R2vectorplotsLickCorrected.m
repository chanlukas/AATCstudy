

clear all   

cd('T:\jan\Collabo Data\HPCpaperPreProcessed')
load AATCR2_Sua_Psth_1ms
cd('T:\jan\Collabo Data\PFCpaperPreProcessed')
load('LickEvokedIndxR2.mat')
% AATCR2_Sua_Psth(:,find(ismember(AnimalCounter,[1,3])),[1,2,3])=AATCR2_Sua_Psth(:,find(ismember(AnimalCounter,[1,3])),[3,2,1]);

AATCR2_Sua_Psth=AATCR2_Sua_Psth(1:5000,:,:);
Bin=25;
AATCR2_Sua_Psth=squeeze(nanmean(reshape(AATCR2_Sua_Psth,Bin,size(AATCR2_Sua_Psth,1)/Bin,size(AATCR2_Sua_Psth,2),size(AATCR2_Sua_Psth,3))));
bc_Psths=AATCR2_Sua_Psth-mean(AATCR2_Sua_Psth(1:1000/Bin,:,:));

for i=1:size(bc_Psths,2)
   bc_Psths(:,i,:) = smoothdata(bc_Psths(:,i,:),'gaussian',25);
end

counter=1;
scounter=1;
for a=1:size(unique(AnimalCounter),2)
for s=1:size(unique(SessionCounter(AnimalCounter==a)),2);

PsthTrace=bc_Psths(:,AnimalCounter==a&SessionCounter==s&Shank==1&LickUp==0&LickDown==0,:);
cellnum=size(find(SessionCounter==s&AnimalCounter==a&Shank==1&LickUp==0&LickDown==0),2);

FRdiff(counter,:)=mean(PsthTrace(:,:,1),2)-mean(PsthTrace(:,:,3),2);
    %for every Bin, calculate distance between population vector in
    %Rewarded vs Unrewarded Contions

for i=1:size(AATCR2_Sua_Psth,1)
    in=squeeze(PsthTrace(i,:,[1,2]))';
out = squareform(pdist(in));
Distance12(i,counter)=out(1,2)./cellnum;
end
for i=1:size(AATCR2_Sua_Psth,1)
    in=squeeze(PsthTrace(i,:,[1,3]))';
out = squareform(pdist(in));
Distance13(i,counter)=out(1,2)./cellnum;
end
for i=1:size(AATCR2_Sua_Psth,1)
    in=squeeze(PsthTrace(i,:,[2,3]))';
out = squareform(pdist(in));
Distance23(i,counter)=out(1,2)./cellnum;
end

if max(LearnedCounter(AnimalCounter==a&SessionCounter==s&TaskModCounter==1))==1&max(Learned2Counter(AnimalCounter==a&SessionCounter==s&TaskModCounter==1))==1;
Late(counter)=1;
else
Late(counter)=0;
end
counter=counter+1
end
end


%%
time=-1+Bin/1000:Bin/1000:4;

    DistanceBC12=(Distance12-mean(Distance12(1:(1000/Bin)-1,:)))./nanstd(Distance12(1:(1000/Bin)-1,:));
    DistanceBC13=(Distance13-mean(Distance13(1:(1000/Bin)-1,:)))./nanstd(Distance13(1:(1000/Bin)-1,:));
    DistanceBC23=(Distance23-mean(Distance23(1:(1000/Bin)-1,:)))./nanstd(Distance23(1:(1000/Bin)-1,:));
    
    
    an12=cat(2,mean(DistanceBC12(1:1000/Bin,Late==1))',mean(DistanceBC12(3000/Bin:4000/Bin,Late==1))');
    an13=cat(2,mean(DistanceBC13(1:1000/Bin,Late==1))',mean(DistanceBC13(3000/Bin:4000/Bin,Late==1))');
    an23=cat(2,mean(DistanceBC23(1:1000/Bin,Late==1))',mean(DistanceBC23(3000/Bin:4000/Bin,Late==1))');
    an=cat(1,an12,an13,an23);

%     [p,tbl,stats] = anova2(an,size(an12,1))

[h,p]=ttest2(an13(:,2),an12(:,2))
[h,p]=ttest2(an13(:,2),an23(:,2))
[h,p]=ttest2(an12(:,2),an23(:,2))


figure()
fig = gcf
fig.Renderer='Painters';
subplot(1,4,[1,2,3])

    DistanceBC=(Distance13-nanmean(Distance13(1:(1000/Bin)-1,:)))./nanstd(Distance13(1:(1000/Bin)-1,:));;
     pp(1,:)=nanmean(DistanceBC(:,Late==1),2);
    errBar1=repmat((std(DistanceBC(:,Late==1)')/sqrt(size(find(Late==1),2))),2,1);
    shadePlot2(time,pp(1,:),errBar1,[1 0 0],.5,'-',3)
%     DistanceBC=(Distance12-nanmean(Distance12(1:(1000/Bin)-1,:)))./nanstd(Distance12(1:(1000/Bin)-1,:));;
%      pp(1,:)=nanmean(DistanceBC(:,Late==1),2);
%     R1R2PFC=DistanceBC;
%     errBar1=repmat((std(DistanceBC(:,Late==1)')/sqrt(size(find(Late==1),2))),2,1);
%     shadePlot2(time,pp(1,:),errBar1,[0 .7 1],.5,'-',3)
%   DistanceBC=(Distance23-nanmean(Distance23(1:(1000/Bin)-1,:)))./nanstd(Distance23(1:(1000/Bin)-1,:));;
%      pp(1,:)=nanmean(DistanceBC(:,Late==1),2);
%     R1R2PFC=DistanceBC;
%     errBar1=repmat((std(DistanceBC(:,Late==1)')/sqrt(size(find(Late==1),2))),2,1);
%     shadePlot2(time,pp(1,:),errBar1,[0 .5 1],.5,'-',3)
    hold on
    plot([0 0],[min(min(pp))-min(min(errBar1)) max(max(pp))+max(max(errBar1))],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
plot([2 2],[min(min(pp))-min(min(errBar1)) max(max(pp))+max(max(errBar1))],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
plot([3 3],[min(min(pp))-min(min(errBar1)) max(max(pp))+max(max(errBar1))],'LineWidth',2,'LineStyle',':','Color',[0 0 0])
box off
axis tight
xlabel('Time [s]')
ylabel('Distance Units')
set(gca,'LineWidth',5)
set(gca,'FontSize',25)

subplot(1,4,[4])


data{1}=an13(:,2);


y=[nanmean(data{1});];

e=[nanstd(data{1})/sqrt(length(data{1}))]

[p1,h]=signrank(data{1})

p=[p1]
c=[0 0 1];
superbar(y,'E',e,'P',p,'BarFaceColor',c)
ylabel('Distance/Cell')
set(gca,'FontSize',20);
set(gca,'LineWidth',5);

[p,h]=ranksum(mean(DistanceBC13(1:1000/Bin,Late==1)),mean(DistanceBC13(3000/Bin:4000/Bin,Late==1)))

%%
clear all   

cd('T:\jan\Collabo Data\HPCpaperPreProcessed')
load AATCR2_Sua_Psth_1ms

cd('T:\jan\Collabo Data\PFCpaperPreProcessed')
load('LickEvokedIndxR2.mat')

AATCR2_Sua_Psth(:,find(ismember(AnimalCounter,[1,3])),[1,2,3])=AATCR2_Sua_Psth(:,find(ismember(AnimalCounter,[1,3])),[3,2,1]);

AATCR2_Sua_Psth=AATCR2_Sua_Psth(1:5000,:,:);

Bin=25;

AATCR2_Sua_Psth=squeeze(sum(reshape(AATCR2_Sua_Psth,Bin,size(AATCR2_Sua_Psth,1)/Bin,size(AATCR2_Sua_Psth,2),size(AATCR2_Sua_Psth,3))));
bc_Psths=AATCR2_Sua_Psth-mean(AATCR2_Sua_Psth(1:1000/Bin,:,:));
norm_bc_Psths=bc_Psths./max(bc_Psths);


for i=1:size(bc_Psths,2)
   bc_Psths(:,i,:) = smoothdata(bc_Psths(:,i,:),'gaussian',25);
end


counter=1;
scounter=1;
for a=1:size(unique(AnimalCounter),2)
for s=1:size(unique(SessionCounter(AnimalCounter==a)),2);

PsthTrace=bc_Psths(:,AnimalCounter==a&SessionCounter==s&Shank==2&LickUp==0&LickDown==0,:);

cellnum=size(find(SessionCounter==s&AnimalCounter==a&Shank==2&LickUp==0&LickDown==0),2);
FRdiff(counter,:)=mean(PsthTrace(:,:,1),2)-mean(PsthTrace(:,:,3),2);

    %for every Bin, calculate distance between population vector in
    %Rewarded vs Unrewarded Contions

for i=1:size(AATCR2_Sua_Psth,1)
    in=squeeze(PsthTrace(i,:,[1,2]))';
out = squareform(pdist(in));
Distance12(i,counter)=out(1,2)./cellnum;
end
for i=1:size(AATCR2_Sua_Psth,1)
    in=squeeze(PsthTrace(i,:,[1,3]))';
out = squareform(pdist(in));
Distance13(i,counter)=out(1,2)./cellnum;
end
for i=1:size(AATCR2_Sua_Psth,1)
    in=squeeze(PsthTrace(i,:,[2,3]))';
out = squareform(pdist(in));
Distance23(i,counter)=out(1,2)./cellnum;
end

if max(LearnedCounter(AnimalCounter==a&SessionCounter==s&TaskModCounter==1))==1&max(Learned2Counter(AnimalCounter==a&SessionCounter==s&TaskModCounter==1))==1;
Late(counter)=1;
else
Late(counter)=0;
end

counter=counter+1
end
end

%%
time=-1+Bin/1000:Bin/1000:4;

    DistanceBC12=(Distance12-mean(Distance12(1:(1000/Bin)-1,:)))./nanstd(Distance12(1:(1000/Bin)-1,:));
    DistanceBC13=(Distance13-mean(Distance13(1:(1000/Bin)-1,:)))./nanstd(Distance13(1:(1000/Bin)-1,:));
    DistanceBC23=(Distance23-mean(Distance23(1:(1000/Bin)-1,:)))./nanstd(Distance23(1:(1000/Bin)-1,:));
    
    
    
    an12=cat(2,mean(DistanceBC12(1:1000/Bin,Late==1))',mean(DistanceBC12(3000/Bin:4000/Bin,Late==1))');
    an13=cat(2,mean(DistanceBC13(1:1000/Bin,Late==1))',mean(DistanceBC13(3000/Bin:4000/Bin,Late==1))');
    an23=cat(2,mean(DistanceBC23(1:1000/Bin,Late==1))',mean(DistanceBC23(3000/Bin:4000/Bin,Late==1))');
    an=cat(1,an12,an13,an23);

%     [p,tbl,stats] = anova2(an,size(an12,1))

[h,p]=ttest2(an13(:,2),an12(:,2))
[h,p]=ttest2(an13(:,2),an23(:,2))
[h,p]=ttest2(an12(:,2),an23(:,2))


figure()
fig = gcf
fig.Renderer='Painters';
subplot(1,4,[1,2,3])

    DistanceBC=(Distance13-nanmean(Distance13(1:(1000/Bin)-1,:)))./nanstd(Distance13(1:(1000/Bin)-1,:));;
     pp(1,:)=nanmean(DistanceBC(:,Late==1),2);
    errBar1=repmat((std(DistanceBC(:,Late==1)')/sqrt(size(find(Late==1),2))),2,1);
    shadePlot2(time,pp(1,:),errBar1,[1 0 0],.5,'-',3)
%     DistanceBC=(Distance12-nanmean(Distance12(1:(1000/Bin)-1,:)))./nanstd(Distance12(1:(1000/Bin)-1,:));;
%      pp(1,:)=nanmean(DistanceBC(:,Late==1),2);
%     R1R2PFC=DistanceBC;
%     errBar1=repmat((std(DistanceBC(:,Late==1)')/sqrt(size(find(Late==1),2))),2,1);
%     shadePlot2(time,pp(1,:),errBar1,[0 .7 1],.5,'-',3)
%   DistanceBC=(Distance23-nanmean(Distance23(1:(1000/Bin)-1,:)))./nanstd(Distance23(1:(1000/Bin)-1,:));;
%      pp(1,:)=nanmean(DistanceBC(:,Late==1),2);
%     R1R2PFC=DistanceBC;
%     errBar1=repmat((std(DistanceBC(:,Late==1)')/sqrt(size(find(Late==1),2))),2,1);
%     shadePlot2(time,pp(1,:),errBar1,[0 .5 1],.5,'-',3)
%     hold on
    plot([0 0],[min(min(pp))-min(min(errBar1)) max(max(pp))+max(max(errBar1))],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
plot([2 2],[min(min(pp))-min(min(errBar1)) max(max(pp))+max(max(errBar1))],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
plot([3 3],[min(min(pp))-min(min(errBar1)) max(max(pp))+max(max(errBar1))],'LineWidth',2,'LineStyle',':','Color',[0 0 0])
box off
axis tight
xlabel('Time [s]')
ylabel('Distance Units')
set(gca,'LineWidth',5)
set(gca,'FontSize',25)


subplot(1,4,[4])


data{1}=an13(:,2);


y=[nanmean(data{1});];

e=[nanstd(data{1})/sqrt(length(data{1}))]

[p1,h]=signrank(data{1})

p=[p1]
c=[0 0 1];
superbar(y,'E',e,'P',p,'BarFaceColor',c)
ylabel('Distance/Cell')
set(gca,'FontSize',20);
set(gca,'LineWidth',5);


