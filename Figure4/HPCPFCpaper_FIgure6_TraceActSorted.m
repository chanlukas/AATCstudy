
clear all   

cd('T:\jan\HeadFixed Data\AATCR2_HPC128-PFC128\AATCR2_Processed_Data')
load AATCR2_Sua_Psth_1ms


AATC_Sua_Psth=AATCR2_Sua_Psth(:,:,:);


%% Pre Process Single Cell Psths

Bin=25;
baseline=1:1000/Bin;

AATC_Sua_PsthBined=squeeze(sum(reshape(AATC_Sua_Psth,Bin,size(AATC_Sua_Psth,1)/Bin,size(AATC_Sua_Psth,2),size(AATC_Sua_Psth,3))));
norm_bc_Psths=(AATC_Sua_PsthBined-nanmean(AATC_Sua_PsthBined(baseline,:,:)))./nanstd(AATC_Sua_PsthBined(baseline,:,:));  %in herz

MinAno=min(AnovaCounter,[],2);

Condition=find(LearnedCounter==1&Learned2Counter==1&ShankS==1);

Psths1=norm_bc_Psths(:,Condition,:);

for i=1:length(Condition)
   Psths(:,i,:) = smoothdata(Psths1(:,i,:),'gaussian',25);
end

idxMD=mean(norm_bc_Psths(3000/Bin:4000/Bin,Condition,1));
[q,sIdxMD]=sort(idxMD);

subplot(1,3,1)
title({'CS+A Trials'})
hold on
imagesc(Psths(:,sIdxMD,1)',[-1 1])
xticklabels = -1:1:4;
xticks = linspace(1, size(Psths, 1), numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
hold on
plot([1000/Bin*1 1000/Bin*1],[size(Condition,2) 0],'LineWidth',2,'LineStyle','--','Color',[1 1 1])
plot([1000/Bin*3 1000/Bin*3],[size(Condition,2) 0],'LineWidth',2,'LineStyle','--','Color',[1 1 1])
plot([1000/Bin*4 1000/Bin*4],[size(Condition,2) 0],'LineWidth',2,'LineStyle',':','Color',[0 0 0])
box off
ylabel('Cell #')
axis tight
set(gca,'LineWidth',3)
set(gca,'FontSize',25)

idxMD=mean(norm_bc_Psths(3000/Bin:4000/Bin,Condition,3));
[q,sIdxMD]=sort(idxMD);


subplot(1,3,2)
title({'CS+B Trials'})
hold on
imagesc(Psths(:,sIdxMD,3)',[-1 1])
set(gca,'YDir','normal')
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
hold on
plot([1000/Bin*1 1000/Bin*1],[size(Condition,2) 0],'LineWidth',2,'LineStyle','--','Color',[1 1 1])
plot([1000/Bin*3 1000/Bin*3],[size(Condition,2) 0],'LineWidth',2,'LineStyle','--','Color',[1 1 1])
plot([1000/Bin*4 1000/Bin*4],[size(Condition,2) 0],'LineWidth',2,'LineStyle',':','Color',[0 0 0])
box off
ylabel('Cell #')
axis tight
set(gca,'LineWidth',3)
set(gca,'FontSize',25)

idxMD=mean(norm_bc_Psths(3000/Bin:4000/Bin,Condition,2));
[q,sIdxMD]=sort(idxMD);

subplot(1,3,3)
title({'CS- Trials'})
hold on
imagesc(Psths(:,sIdxMD,2)',[-1 1])
set(gca,'YDir','normal')
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
hold on
plot([1000/Bin*1 1000/Bin*1],[size(Condition,2) 0],'LineWidth',2,'LineStyle','--','Color',[1 1 1])
plot([1000/Bin*3 1000/Bin*3],[size(Condition,2) 0],'LineWidth',2,'LineStyle','--','Color',[1 1 1])
box off
ylabel('Cell #')
axis tight
set(gca,'LineWidth',3)
set(gca,'FontSize',25)

clear all   

cd('T:\jan\HeadFixed Data\AATCR2_HPC128-PFC128\AATCR2_Processed_Data')
load AATCR2_Sua_Psth_1ms


AATC_Sua_Psth=AATCR2_Sua_Psth(:,:,:);


%% Pre Process Single Cell Psths

Bin=25;
baseline=1:1000/Bin;

AATC_Sua_PsthBined=squeeze(sum(reshape(AATC_Sua_Psth,Bin,size(AATC_Sua_Psth,1)/Bin,size(AATC_Sua_Psth,2),size(AATC_Sua_Psth,3))));
norm_bc_Psths=(AATC_Sua_PsthBined-nanmean(AATC_Sua_PsthBined(baseline,:,:)))./nanstd(AATC_Sua_PsthBined(baseline,:,:));  %in herz

MinAno=min(AnovaCounter,[],2);

Condition=find(LearnedCounter==1&Learned2Counter==1&ShankS==2);

Psths1=norm_bc_Psths(:,Condition,:);

for i=1:length(Condition)
   Psths(:,i,:) = smoothdata(Psths1(:,i,:),'gaussian',25);
end

idxMD=mean(norm_bc_Psths(3000/Bin:4000/Bin,Condition,1));
[q,sIdxMD]=sort(idxMD);

figure()
subplot(1,3,1)
title({'CS+A Trials'})
hold on
imagesc(Psths(:,sIdxMD,1)',[-1 1])
xticklabels = -1:1:4;
xticks = linspace(1, size(Psths, 1), numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
hold on
plot([1000/Bin*1 1000/Bin*1],[size(Condition,2) 0],'LineWidth',2,'LineStyle','--','Color',[1 1 1])
plot([1000/Bin*3 1000/Bin*3],[size(Condition,2) 0],'LineWidth',2,'LineStyle','--','Color',[1 1 1])
plot([1000/Bin*4 1000/Bin*4],[size(Condition,2) 0],'LineWidth',2,'LineStyle',':','Color',[0 0 0])
box off
ylabel('Cell #')
axis tight
set(gca,'LineWidth',3)
set(gca,'FontSize',25)

idxMD=mean(norm_bc_Psths(3000/Bin:4000/Bin,Condition,3));
[q,sIdxMD]=sort(idxMD);


subplot(1,3,2)
title({'CS+B Trials'})
hold on
imagesc(Psths(:,sIdxMD,3)',[-1 1])
set(gca,'YDir','normal')
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
hold on
plot([1000/Bin*1 1000/Bin*1],[size(Condition,2) 0],'LineWidth',2,'LineStyle','--','Color',[1 1 1])
plot([1000/Bin*3 1000/Bin*3],[size(Condition,2) 0],'LineWidth',2,'LineStyle','--','Color',[1 1 1])
plot([1000/Bin*4 1000/Bin*4],[size(Condition,2) 0],'LineWidth',2,'LineStyle',':','Color',[0 0 0])
box off
ylabel('Cell #')
axis tight
set(gca,'LineWidth',3)
set(gca,'FontSize',25)

idxMD=mean(norm_bc_Psths(3000/Bin:4000/Bin,Condition,2));
[q,sIdxMD]=sort(idxMD);

subplot(1,3,3)
title({'CS- Trials'})
hold on
imagesc(Psths(:,sIdxMD,2)',[-1 1])
set(gca,'YDir','normal')
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
hold on
plot([1000/Bin*1 1000/Bin*1],[size(Condition,2) 0],'LineWidth',2,'LineStyle','--','Color',[1 1 1])
plot([1000/Bin*3 1000/Bin*3],[size(Condition,2) 0],'LineWidth',2,'LineStyle','--','Color',[1 1 1])
box off
ylabel('Cell #')
axis tight
set(gca,'LineWidth',3)
set(gca,'FontSize',25)
