clear all
cd('/Users/Jan/Documents/PFC Paper/AATCdata/PreprocessedDATA_PFCpaper')

load('HPCPFC_1thLicksInNoutOfTrials.mat')
load('AATC_PFC_Sua_Psth_1ms.mat')
% load('TraceCounters.mat')
AATC_Sua_Psth=AATC_PFC_Sua_Psth;

Bin=25;
baseline=1:1000/Bin;
AATC_Sua_Psth=LickTrigPsthInsideTrialPFC';

AATC_Sua_PsthBined=squeeze(mean(reshape(AATC_Sua_Psth,Bin,size(AATC_Sua_Psth,1)/Bin,size(AATC_Sua_Psth,2))));

PreTrialBaseline=nanmean(AATC_Sua_PsthBined(baseline,:,1));
PreTrialBaselineSTD=nanstd(AATC_Sua_PsthBined(baseline,:,1));
PreTrialBaselineSTD(find(PreTrialBaselineSTD==0))=NaN;


bc_Psths=(AATC_Sua_PsthBined-PreTrialBaseline)./PreTrialBaselineSTD; %in herz
for i=1:size(bc_Psths,2)
   norm_bc_Psths1(:,i,:) = smoothdata(bc_Psths(:,i,:),'gaussian',20);
end

Condition=find(isnan(mean(norm_bc_Psths1))==0&TrgDayCounter<3&LearnedCounter==0);
Condition2=find(isnan(mean(norm_bc_Psths1))==0&LearnedCounter==1);

figure()


subplot(1,2,1)

[v,i]=sort(mean(norm_bc_Psths1(1750/Bin:2250/Bin,Condition2)))

imagesc(norm_bc_Psths1(:,Condition2(i))',[-1 1])
set(gca,'YDir','normal')
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
hold on
plot([2000/Bin*1 2000/Bin*1],[size(Condition2,2) 0],'LineWidth',2,'LineStyle','--','Color',[1 1 1])
box off
ylabel('Cell #')
% axis tight
set(gca,'LineWidth',3)
set(gca,'FontSize',25)

Tresh=1;
EvR=squeeze(mean((norm_bc_Psths1(1750/Bin:2250/Bin,:,:))));
LickUpPFC(1:size(AATC_Sua_Psth,2))=0;
LickUpPFC(find(EvR>PreTrialBaselineSTD*Tresh))=1;
Tresh=-1;
LickDownPFC(1:size(AATC_Sua_Psth,2))=0;
LickDownPFC(find(EvR<PreTrialBaselineSTD*Tresh))=1;

subplot(1,2,2)
time=-2+0.001*Bin:.001*Bin:2;

pp(1,:)=nanmean(norm_bc_Psths1(:,Condition2)');
errBar=repmat(nanstd(norm_bc_Psths1(:,Condition2)')/sqrt(size(Condition2,2)),2,1);
shadePlot2(time,pp(1,:),errBar(:,:),[1 0 0],.3)
hold on

pp(1,:)=nanmean(norm_bc_Psths1(:,isnan(mean(norm_bc_Psths1))==0&LearnedCounter==1&LickDownPFC==1)');
errBar=repmat(nanstd(norm_bc_Psths1(:,isnan(mean(norm_bc_Psths1))==0&LearnedCounter==1&LickDownPFC==1)')/sqrt(size(find(isnan(mean(norm_bc_Psths1))==0&LearnedCounter==1&LickDownPFC==1),2)),2,1);
shadePlot2(time,pp(1,:),errBar(:,:),[1 0 0],.3)

pp(1,:)=nanmean(norm_bc_Psths1(:,isnan(mean(norm_bc_Psths1))==0&LearnedCounter==1&LickUpPFC==1)');
errBar=repmat(nanstd(norm_bc_Psths1(:,isnan(mean(norm_bc_Psths1))==0&LearnedCounter==1&LickUpPFC==1)')/sqrt(size(find(isnan(mean(norm_bc_Psths1))==0&LearnedCounter==1&LickUpPFC==1),2)),2,1);
shadePlot2(time,pp(1,:),errBar(:,:),[1 0 0],.3)

plot([0 0],[max(max(pp))+max(max(errBar)) min(min(pp))-min(min(errBar))],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
box off
% axis tight% axis tight


length(find(LickUpPFC==1))
length(find(LickDownPFC==1))

length(find(LickUpPFC==1))/length(find(isnan(mean(norm_bc_Psths1))==0&LearnedCounter==1))

length(find(LickDownPFC==1))



%%
clearvars -except LickUpPFC LickDownPFC
% cd('T:\jan\Collabo Data\PFCpaperPreProcessed')

load('HPCPFC_1thLicksInNoutOfTrials.mat')
% cd('T:\jan\Collabo Data\HPCpaperPreProcessed')

load('AATC_Sua_Psth_1ms.mat')

% cd('T:\jan\Collabo Data\PFCpaperPreProcessed')
% load('TraceCounters.mat')
Bin=25;
baseline=1:1000/Bin;
AATC_Sua_Psth=LickTrigPsthInsideTrialHPC';

AATC_Sua_PsthBined=squeeze(mean(reshape(AATC_Sua_Psth,Bin,size(AATC_Sua_Psth,1)/Bin,size(AATC_Sua_Psth,2))));

PreTrialBaseline=nanmean(AATC_Sua_PsthBined(baseline,:,1));
PreTrialBaselineSTD=nanstd(AATC_Sua_PsthBined(baseline,:,1));
PreTrialBaselineSTD(find(PreTrialBaselineSTD==0))=NaN;


bc_Psths=(AATC_Sua_PsthBined-PreTrialBaseline)./PreTrialBaselineSTD; %in herz
for i=1:size(bc_Psths,2)
   norm_bc_Psths1(:,i,:) = smoothdata(bc_Psths(:,i,:),'gaussian',20);
end

Condition=find(isnan(mean(norm_bc_Psths1))==0&TrgDayCounter<3&LearnedCounter==0);
Condition2=find(isnan(mean(norm_bc_Psths1))==0&LearnedCounter==1);

figure()

subplot(1,2,1)

[v,i]=sort(mean(norm_bc_Psths1(1750/Bin:2250/Bin,Condition2)))

imagesc(norm_bc_Psths1(:,Condition2(i))',[-1 1])
set(gca,'YDir','normal')
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
hold on
plot([2000/Bin*1 2000/Bin*1],[size(Condition2,2) 0],'LineWidth',2,'LineStyle','--','Color',[1 1 1])
box off
ylabel('Cell #')
% axis tight
set(gca,'LineWidth',3)
set(gca,'FontSize',25)


subplot(1,2,2)
time=-2+0.001*Bin:.001*Bin:2;

baseline=1:1000/Bin;
PreTrialBaseline=nanmean(AATC_Sua_PsthBined(baseline,:,1));
PreTrialBaselineSTD=nanstd(AATC_Sua_PsthBined(baseline,:,1));
Tresh=1;
EvR=squeeze(mean((norm_bc_Psths1(1750/Bin:2250/Bin,:,:))));
LickUpHPC(1:size(AATC_Sua_Psth,2))=0;
LickUpHPC(find(EvR>PreTrialBaselineSTD*Tresh))=1;
Tresh=-1;
LickDownHPC(1:size(AATC_Sua_Psth,2))=0;
LickDownHPC(find(EvR<PreTrialBaselineSTD*Tresh))=1;

pp(1,:)=nanmean(norm_bc_Psths1(:,Condition2)');
errBar=repmat(std(norm_bc_Psths1(:,Condition2)')/sqrt(size(Condition2,2)),2,1);
shadePlot2(time,pp(1,:),errBar(:,:),[1 0 0],.3)
hold on

pp(1,:)=nanmean(norm_bc_Psths1(:,isnan(mean(norm_bc_Psths1))==0&LearnedCounter==1&LickDownHPC==1)');
errBar=repmat(std(norm_bc_Psths1(:,isnan(mean(norm_bc_Psths1))==0&LearnedCounter==1&LickDownHPC==1)')/sqrt(size(find(isnan(mean(norm_bc_Psths1))==0&LearnedCounter==1&LickDownHPC==1),2)),2,1);
shadePlot2(time,pp(1,:),errBar(:,:),[1 0 0],.3)

pp(1,:)=nanmean(norm_bc_Psths1(:,isnan(mean(norm_bc_Psths1))==0&LearnedCounter==1&LickUpHPC==1)');
errBar=repmat(std(norm_bc_Psths1(:,isnan(mean(norm_bc_Psths1))==0&LearnedCounter==1&LickUpHPC==1)')/sqrt(size(find(isnan(mean(norm_bc_Psths1))==0&LearnedCounter==1&LickUpHPC==1),2)),2,1);
shadePlot2(time,pp(1,:),errBar(:,:),[1 0 0],.3)

plot([0 0],[max(max(pp))+max(max(errBar)) min(min(pp))-min(min(errBar))],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
box off
% axis tight

length(find(LickUpHPC==1&isnan(mean(norm_bc_Psths))==0&LearnedCounter==1))/length(find(isnan(mean(norm_bc_Psths))==0&LearnedCounter==1))
length(find(LickDownHPC==1&isnan(mean(norm_bc_Psths))==0&LearnedCounter==1))/length(find(isnan(mean(norm_bc_Psths))==0&LearnedCounter==1))


length(find(LickUpHPC==1))
length(find(LickDownHPC==1))
% save ('LickEvokedIndx','LickUpHPC','LickDownHPC','LickUpPFC', 'LickDownPFC')
% 
