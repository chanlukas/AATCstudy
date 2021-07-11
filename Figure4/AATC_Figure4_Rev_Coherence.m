clear all

cd('/Users/Jan/Documents/PFC Paper/AATCdata/PreprocessedDATA_PFCpaper')



load('HPCPFC_coherence.mat')
load('HPCPFC_TFA.mat')

Pow=PowPFC;

Hz50=zeros(1,21);

Hz50(find(squeeze(nanmean(Pow(50,1:5999,:,2),2))>150))=1;

condition=find(Learned==1&Hz50==0)
condition2=find(Learned==0&Hz50==0)
 
tt=1:131 
figure();
fig = gcf
fig.Renderer='Painters';
subplot(1,3,[1])
imagesc(tt,f,squeeze(nanmean(CoherenceR(tt,:,condition),3))',[.5 .8])
hold on
yticklabels = 0:30:150;
yticks = linspace(1, size([1:150],2), numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
plot([31 31],[1 150],'LineWidth',2,'LineStyle',':','Color',[1 1 1])
plot([51 51],[1 150],'LineWidth',2,'LineStyle',':','Color',[1 1 1])
plot([61 61],[1 150],'LineWidth',2,'LineStyle',':','Color',[1 1 1])
xticklabels = -5:5:25;
xticks = linspace(1, 281, numel(xticklabels));
% set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
box off
ylabel('Frequency')
title('CS+ Trials Avg')
set(gca,'FontSize',25);
set(gca,'LineWidth',5);
set(gca,'YDir','normal')


subplot(1,3,[2])
imagesc(tt,f,squeeze(nanmean(CoherenceU(tt,:,condition),3))',[.5 .8])
yticklabels = 0:30:150;
yticks = linspace(1, size([1:150],2), numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
hold on
plot([31 31],[1 150],'LineWidth',2,'LineStyle','--','Color',[1 1 1])
plot([51 51],[1 150],'LineWidth',2,'LineStyle','--','Color',[1 1 1])
xticklabels = -5:5:25;
xticks = linspace(2, 281, numel(xticklabels));
% set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
box off
ylabel('Frequency')
title('CS+ Trials Avg')
set(gca,'FontSize',25);
set(gca,'LineWidth',5);
set(gca,'YDir','normal')

% subplot(2,2,[3,4])
% 
% meanCo=squeeze(nanmean(CoherenceR(tt,:,condition),3))';
% meanCoU=squeeze(nanmean(CoherenceU(tt,:,condition),3))';
% pp=mean(meanCo([288:509],:));
% errBar=repmat(std(meanCo([288:509],:)),2,1)./sqrt(length(condition));
% shadePlot2(tt,pp,errBar,[0 0 1],.5)
% hold on
% pp=mean(meanCoU([288:509],:));
% errBar=repmat(std(meanCoU([288:509],:)),2,1)./sqrt(length(condition));
% shadePlot2(tt,pp,errBar,[1 0 0],.5)
% hold on
% plot([33 33],[0 1.5],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
% plot([53 53],[0 1.5],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
% plot([63 63],[0 1.5],'LineWidth',2,'LineStyle',':','Color',[0 0 0])
% xticklabels = -5:5:25;
% xticks = linspace(1, 281, numel(xticklabels));
% % set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
% box off
% axis tight
% xlabel('Time [s]')
% ylabel('Coherence')
% title('80-150 Hz')
% set(gca,'FontSize',25);
% set(gca,'LineWidth',5);
% ylim([0 1])

subplot(1,3,[3])

CoherenceMeanR=squeeze(nanmean(CoherenceR(51:60,:,:),1))';
CoherenceMeanU=squeeze(nanmean(CoherenceU(51:60,:,:),1))';

pp=nanmean(CoherenceMeanR);
errBar=repmat(nanstd(CoherenceMeanR),2,1)./sqrt(length(condition));
shadePlot2(f,pp,errBar,[0 0 1],.5)
hold on
pp=nanmean(CoherenceMeanU);
errBar=repmat(nanstd(CoherenceMeanU),2,1)./sqrt(length(condition));
shadePlot2(f,pp,errBar,[1 0 0],.5)

for i=1:1:509
    CoherenceNan=cat(2,CoherenceMeanR(condition,i),CoherenceMeanU(condition,i));
    x= sum(isnan(CoherenceNan),2);
    CoherenceCleanR=squeeze(CoherenceNan(x==0,1));
    CoherenceCleanU=squeeze(CoherenceNan(x==0,2));
[ppp(i,1)]=permutationTest(CoherenceCleanU',CoherenceCleanR',1000);
clearvars CoherenceNan x CoherenceCleanR CoherenceCleanU
end

idx=(find(ppp<0.05))
x=NaN(1,509)
x(find(ppp<0.05))=.75;
c=plot(f,x,'-k')

% figure()
% data{1}=nanmean(squeeze(nanmean(CoherenceR(12:41,32:509,condition))));
% data{2}=nanmean(squeeze(nanmean(CoherenceU(12:41,32:509,condition))));
% 
% y=[nanmean(data{1});nanmean(data{2})];
% 
% e=[nanstd(data{1})/sqrt(length(data{1})),...
%     nanstd(data{2})/sqrt(length(data{2}))]
% 
% [p1,h]=ranksum(data{1},data{2})
% 
% p=[ NaN p1;
%     p1,NaN ]
% 
% c=[0 0 1;1 0 0]
% superbar(y,'E',e,'P',p,'BarFaceColor',c)
% names={'CS+';'CS-'}
% set(gca,'Xtick',1:1:2)
% set(gca,'xticklabel',names)
% ylabel('')
% title('10-150 Hz')
% set(gca,'FontSize',20);
% set(gca,'LineWidth',5);
% clear all
% 
% cd('T:\jan\HeadFixed Data\AATC\AATC_Processed_Data') 
% 
% 
% load AATC_PFC_cell_HPC_Lfp_PhaseLocking 
% 
% data=cat(1,MRLHighGammaR,MRLHighGammaU,MRLGammaR,MRLGammaU,MRLThetaR,MRLThetaU);
% data=mat2cell(data,6,1016)
% figure()
% 
% data{1}=MRLThetaR(Learned==0);
% data{2}=MRLThetaU(Learned==0);
% data{3}=MRLThetaR(Learned==1);
% data{4}=MRLThetaU(Learned==1);
% 
% data{5}=MRLGammaR(Learned==0);
% data{6}=MRLGammaU(Learned==0);
% data{7}=MRLGammaR(Learned==1);
% data{8}=MRLGammaU(Learned==1);
% 
% data{9}=MRLHighGammaR(Learned==0);
% data{10}=MRLHighGammaU(Learned==0);
% data{11}=MRLHighGammaR(Learned==1);
% data{12}=MRLHighGammaU(Learned==1);
% 
% y=[nanmean(data{1}),nanmean(data{2}),...
%    nanmean(data{3}),nanmean(data{4}),...
%    nanmean(data{5}),nanmean(data{6}),...
%    nanmean(data{7}),nanmean(data{8}),...
%    nanmean(data{9}),nanmean(data{10}),...
%    nanmean(data{11}),nanmean(data{12})]
% 
% e=[nanstd(data{1})/sqrt(length(data{1})),...
%     nanstd(data{2})/sqrt(length(data{2})),...
%     nanstd(data{3})/sqrt(length(data{3})),...
%     nanstd(data{4})/sqrt(length(data{4})),...
%     nanstd(data{5})/sqrt(length(data{5})),...
%     nanstd(data{6})/sqrt(length(data{6})),...
%     nanstd(data{7})/sqrt(length(data{7})),...
%     nanstd(data{8})/sqrt(length(data{8})),...
%     nanstd(data{9})/sqrt(length(data{9})),...
%     nanstd(data{10})/sqrt(length(data{10})),...
%     nanstd(data{11})/sqrt(length(data{11})),...
%     nanstd(data{12})/sqrt(length(data{12}))]
% 
% [h,p1]=ttest2(data{3},data{4});
% [h,p2]=ttest2(data{7},data{8});
% [h,p3]=ttest2(data{11},data{12});
% 
% 
% 
% 
% p=[NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN;
% NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN;
% NaN NaN NaN p1 NaN NaN NaN NaN NaN NaN NaN NaN;
% NaN NaN p1 NaN NaN NaN NaN NaN NaN NaN NaN NaN;
% NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN;
% NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN;
% NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN;
% NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN;
% NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN;
% NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN;
% NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN;
% NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN;]
% 
% c=[0 0 1;1 0 0;0 0 1;1 0 0;0 0 1;1 0 0;0 0 1;1 0 0;0 0 1;1 0 0;0 0 1;1 0 0]
% superbar(y,'E',e,'P',p,'BarFaceColor',c)
% names={'CS+';'CS-'}
% set(gca,'Xtick',1:1:2)
% set(gca,'xticklabel',names)
% ylabel('MRL')
% title('PFC unit-HPC High Gamma')
% 
% set(gca,'FontSize',20);
% set(gca,'LineWidth',5);
% legend boxoff
% 
% clear all
% 
% load AATC_HPC_cell_PFC_Lfp_PhaseLocking 
% 
% % idxHz50=find(Hz50==0)
% % SessionClean(1:length(Session))=0
% % SessionClean(ismember(SessionAll,idxHz50))=1
% data=cat(1,MRLHighGammaR,MRLHighGammaU,MRLGammaR,MRLGammaU,MRLThetaR,MRLThetaU);
% data=mat2cell(data,6,622)
% figure()
% 
% data{1}=MRLThetaR(Learned==0)%&SessionClean==1);
% data{2}=MRLThetaU(Learned==0)%&SessionClean==1);
% data{3}=MRLThetaR(Learned==1)%&SessionClean==1);
% data{4}=MRLThetaU(Learned==1)%&SessionClean==1);
% 
% data{5}=MRLGammaR(Learned==0)%&SessionClean==1);
% data{6}=MRLGammaU(Learned==0)%&SessionClean==1);
% data{7}=MRLGammaR(Learned==1)%&SessionClean==1);
% data{8}=MRLGammaU(Learned==1)%&SessionClean==1);
% 
% data{9}=MRLHighGammaR(Learned==0)%&SessionClean==1);
% data{10}=MRLHighGammaU(Learned==0)%&SessionClean==1);
% data{11}=MRLHighGammaR(Learned==1)%&SessionClean==1);
% data{12}=MRLHighGammaU(Learned==1)%&SessionClean==1);
% 
% y=[nanmean(data{1}),nanmean(data{2}),...
%    nanmean(data{3}),nanmean(data{4}),...
%    nanmean(data{5}),nanmean(data{6}),...
%    nanmean(data{7}),nanmean(data{8}),...
%    nanmean(data{9}),nanmean(data{10}),...
%    nanmean(data{11}),nanmean(data{12})]
% 
% e=[nanstd(data{1})/sqrt(length(data{1})),...
%     nanstd(data{2})/sqrt(length(data{2})),...
%     nanstd(data{3})/sqrt(length(data{3})),...
%     nanstd(data{4})/sqrt(length(data{4})),...
%     nanstd(data{5})/sqrt(length(data{5})),...
%     nanstd(data{6})/sqrt(length(data{6})),...
%     nanstd(data{7})/sqrt(length(data{7})),...
%     nanstd(data{8})/sqrt(length(data{8})),...
%     nanstd(data{9})/sqrt(length(data{9})),...
%     nanstd(data{10})/sqrt(length(data{10})),...
%     nanstd(data{11})/sqrt(length(data{11})),...
%     nanstd(data{12})/sqrt(length(data{12}))]
% 
% [h,p1]=ttest2(data{3},data{4})
% [h,p2]=ttest2(data{7},data{8})
% [h,p3]=ttest2(data{11},data{12})
% 
% 
% 
% 
% p=[NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN;
% NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN;
% NaN NaN NaN p1 NaN NaN NaN NaN NaN NaN NaN NaN;
% NaN NaN p1 NaN NaN NaN NaN NaN NaN NaN NaN NaN;
% NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN;
% NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN;
% NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN;
% NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN;
% NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN;
% NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN;
% NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN;
% NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN;]
% 
% c=[0 0 1;1 0 0;0 0 1;1 0 0;0 0 1;1 0 0;0 0 1;1 0 0;0 0 1;1 0 0;0 0 1;1 0 0]
% superbar(y,'E',e,'P',p,'BarFaceColor',c)
% names={'CS+';'CS-'}
% set(gca,'Xtick',1:1:2)
% set(gca,'xticklabel',names)
% ylabel('MRL')
% title('HPC unit-PFC Theta')
% 
% set(gca,'FontSize',20);
% set(gca,'LineWidth',5);
% legend boxoff
% 
% 
% 
% 
