
clear all      
cd('T:\jan\Collabo Data\HPCpaperPreProcessed')
load AATC_Sua_Psth_1ms




%% Pre Process Single Cell Psths

Bin=50;
baseline=1:1000/Bin;

AATC_Sua_PsthBined=squeeze(sum(reshape(AATC_Sua_Psth,Bin,size(AATC_Sua_Psth,1)/Bin,size(AATC_Sua_Psth,2),size(AATC_Sua_Psth,3))));
norm_bc_Psths=(AATC_Sua_PsthBined-nanmean(AATC_Sua_PsthBined(baseline,:,:)))./nanstd(AATC_Sua_PsthBined(baseline,:,:));  %in herz

MinAno=min(AnovaCounter,[],2);

Condition=find(TrgDayCounter<3&LearnedCounter==0);

Psths1=norm_bc_Psths(:,Condition,:);

for i=1:length(Condition)
   Psths(:,i,:) = smoothdata(Psths1(:,i,:),'gaussian',5);
end


idxMD=mean(norm_bc_Psths(3000/Bin:4000/Bin,Condition,1));
[q,sIdxMD]=sort(idxMD);

figure()
% subplot(1,4,2)
subplot(1,2,1)

title({'Post Learning CS+ Trials'})
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

idxMD=mean(norm_bc_Psths(3000/Bin:4000/Bin,Condition,2));
[q,sIdxMD]=sort(idxMD);


% subplot(1,4,3)
subplot(1,2,2)

title({'Post Learning CS- Trials'})
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



% % 
% subplot(1,4,4)
% 
% load AATC_Sua_RipPsthNew
% Bin=5;
% baseline=[1:250/Bin];
% AATC_Sua_RipPsthBined=squeeze(nanmean(reshape(AATC_Sua_RipPsth,Bin,size(AATC_Sua_RipPsth,1)/Bin,size(AATC_Sua_RipPsth,2),size(AATC_Sua_RipPsth,3))));
% bc_Psths2=(AATC_Sua_RipPsthBined-nanmean(AATC_Sua_RipPsthBined(baseline,:)));  %in herz
% imagesc(bc_Psths2(400/Bin:600/Bin,sIdxMD,1)',[-5 25])
% set(gca,'YDir','normal')
% set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
% 
%% Analysis sketches
clear all           
cd('T:\jan\Collabo Data\PFCpaperPreProcessed')
load AATC_PFC_Sua_Psth_1ms


AATC_Sua_Psth=AATC_PFC_Sua_Psth;


%% Pre Process Single Cell Psths

Bin=50;
baseline=1:1000/Bin;

AATC_Sua_PsthBined=squeeze(sum(reshape(AATC_Sua_Psth,Bin,size(AATC_Sua_Psth,1)/Bin,size(AATC_Sua_Psth,2),size(AATC_Sua_Psth,3))));
norm_bc_Psths=(AATC_Sua_PsthBined-nanmean(AATC_Sua_PsthBined(baseline,:,:)))./nanstd(AATC_Sua_PsthBined(baseline,:,:));  %in herz

MinAno=min(AnovaCounter,[],2);

Condition=find(TrgDayCounter<3&LearnedCounter==0);

Psths1=norm_bc_Psths(:,Condition,:);

for i=1:length(Condition)
   Psths(:,i,:) = smoothdata(Psths1(:,i,:),'gaussian',5);
end


idxMD=mean(norm_bc_Psths(3000/Bin:4000/Bin,Condition,1));
[q,sIdxMD]=sort(idxMD);

figure()
% subplot(1,4,2)
subplot(1,2,1)

title({'Post Learning CS+ Trials'})
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

idxMD=mean(norm_bc_Psths(3000/Bin:4000/Bin,Condition,2));
[q,sIdxMD]=sort(idxMD);


% subplot(1,4,3)
subplot(1,2,2)

title({'Post Learning CS- Trials'})
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



% % 
% subplot(1,4,4)
% 
% load AATC_Sua_RipPsthNew
% Bin=5;
% baseline=[1:250/Bin];
% AATC_Sua_RipPsthBined=squeeze(nanmean(reshape(AATC_Sua_RipPsth,Bin,size(AATC_Sua_RipPsth,1)/Bin,size(AATC_Sua_RipPsth,2),size(AATC_Sua_RipPsth,3))));
% bc_Psths2=(AATC_Sua_RipPsthBined-nanmean(AATC_Sua_RipPsthBined(baseline,:)));  %in herz
% imagesc(bc_Psths2(400/Bin:600/Bin,sIdxMD,1)',[-5 25])
% set(gca,'YDir','normal')
% set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
% 
