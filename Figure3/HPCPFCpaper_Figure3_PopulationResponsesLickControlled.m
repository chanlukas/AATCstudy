
clear all      
cd('T:\jan\Collabo Data\HPCpaperPreProcessed')
load AATC_Sua_Psth_1ms
cd('T:\jan\Collabo Data\PFCpaperPreProcessed')
load('LickEvokedIndx.mat')

Bin=100;
AATC_Sua_Psth=squeeze(mean(reshape(AATC_Sua_Psth,Bin,size(AATC_Sua_Psth,1)/Bin,size(AATC_Sua_Psth,2),size(AATC_Sua_Psth,3))));
bc_Psths=AATC_Sua_Psth-mean(AATC_Sua_Psth(1:1000/Bin,:,:));


% for i=1:size(bc_Psths,2)
%    bc_Psths(:,i,:) = smoothdata(bc_Psths(:,i,:),'gaussian',25);
% end


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
time=-1+Bin/1000:Bin/1000:4;

[Tr] = pca(PostExample(:,:,1)','NumComponents',3);
[Tur] = pca(PostExample(:,:,2)','NumComponents',3);
HPCexamplePopVectorCSplus=Tr;
HPCexamplePopVectorCSminus=Tur;


for i=1:3
   Tr(:,i) = smoothdata(Tr(:,i),'gaussian',7);
   Tur(:,i) = smoothdata(Tur(:,i),'gaussian',7);

end

HPCexamplePopVectorCSplusSmooth=Tr;
HPCexamplePopVectorCSminusSmooth=Tur;

cd('C:\Users\Jan\Desktop')
save('HPCexamplePopVectors','HPCexamplePopVectorCSplus','HPCexamplePopVectorCSminus','HPCexamplePopVectorCSplusSmooth','HPCexamplePopVectorCSminusSmooth')
plot3(Tr(base,1),Tr(base,2),Tr(base,3),'LineWidth',1,'Color',[0 0 1])
hold on
plot3(Tr(Tone,1),Tr(Tone,2),Tr(Tone,3),'LineWidth',1,'Color',[0 0 1])
hold on
plot3(Tr(Trace,1),Tr(Trace,2),Tr(Trace,3),'LineWidth',1,'Color',[0 0 1])
hold on
plot3(Tr(Reward,1),Tr(Reward,2),Tr(Reward,3),'LineWidth',1,'Color',[0 0 1])
hold on

plot3(Tur(base,1),Tur(base,2),Tur(base,3),'LineWidth',1,'Color',[1 0 0])
hold on
plot3(Tur(Tone,1),Tur(Tone,2),Tur(Tone,3),'LineWidth',1,'Color',[1 0 0])
hold on
plot3(Tur(Trace,1),Tur(Trace,2),Tur(Trace,3),'LineWidth',1,'Color',[1 0 0])
hold on
plot3(Tur(Reward,1),Tur(Reward,2),Tur(Reward,3),'LineWidth',1,'Color',[1 0 0])
hold on
grid on
axis tight
hold on

figure()
hold on
    DistanceBC=(Distance-nanmean(Distance(1:(1000/Bin)-1,:)))./nanstd(Distance(1:(1000/Bin)-1,:));
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

[p,h]=signrank(nanmean(DistanceBC(30:40,Late==1)))

cd('T:\jan\Collabo Data\HPCpaperPreProcessed')
load SVM

figure()
n1{1}=1-allMcRb(LearnedCounter==1)
n1{2}=1-allMcRT(LearnedCounter==1)
names={'Baseline','Trace'};
% ttl={'SVM Classifier'};
ttl={''};
colors=[0 1 0;0 0 0;];
ylabel={'Percent Correct'}
style=([1 1 1;1 1 1]);
BeehivePlotRS(n1,names,ylabel,ttl,colors,style)
hold on
plot([.5 2.5],[0.5 0.5],':r','LineWidth',2)
set(gca,'FontSize',25);
set(gca,'LineWidth',5);


mean(n1{2})
[p,h]=signrank(allMcRT(LearnedCounter==1))

clearvars -except DistanceBC Late incl AvgLickCounter

  cd('T:\jan\Collabo Data\HPCpaperPreProcessed')
    load Speed
    [ConFac]=ballCali()
  
  Bin=100*30;

  SpeedBinned=squeeze(nanmean(reshape(Speed,Bin,size(Speed,1)/Bin,2,size(Speed,3))))*(30000/Bin); 
  Speed=SpeedBinned(91:140,:,:);
  SpeedDiff=squeeze(Speed(:,1,:)-Speed(:,2,:));
  SpeedDiffZ=(SpeedDiff-mean(SpeedDiff(1:10,:)))./std(SpeedDiff(1:10,:));

LickPsth=AvgLickCounter(incl==1,:,1);
UnLickPsth=AvgLickCounter(incl==1,:,2);
Bin=100;
LickPsth=permute(LickPsth(:,:),[2,1,3]);
LickBinned=squeeze(sum(reshape(LickPsth,Bin,size(LickPsth,1)/Bin,size(LickPsth,2)))); 
UnLickPsth=permute(UnLickPsth(:,:),[2,1,3]);
UnLickBinned=squeeze(sum(reshape(UnLickPsth,Bin,size(UnLickPsth,1)/Bin,size(UnLickPsth,2)))); 
  
LickPsthInc=LickBinned(:,:);
UnLickPsthInc=UnLickBinned(:,:);
x=(LickPsthInc-UnLickPsthInc)
x=(x-nanmean(x(1:10,:)))./nanstd(x(1:10,:));
y=(DistanceBC(:,:)-(mean(DistanceBC(1:10,:))))./std(DistanceBC(1:10,:));
z=SpeedDiffZ(:,incl==1);
% 
% for i=1:length(find(incl==1))
% [h,p]=corr(x(10:40,i),y(10:40,i))
% Pcounter(i)=p;
% end
% SigCorPopLick=length(find(Pcounter<.05))

x(x==inf|x==-inf)=1;
x(isnan(x)==1)=1;
[h,p]=corrcoef(nanmean(x(10:40,Late==1)),nanmean(y(10:40,Late==1)))

% for i=1:length(find(incl==1))
% [h,p]=corr(x(10:40,i),z(10:40,i))
% Pcounter(i)=p;
% end
% SigCorPopSpeed=length(find(Pcounter<.05))

z(z==inf)=0
z(isnan(z)==1)=0
[h,p]=corrcoef(nanmean(z(10:40,Late==1)),nanmean(y(10:40,Late==1)))

[v,idx]=sort(Late);
subplot(1,3,1)
imagesc(x(:,Late==1)',[0 25])
box off
subplot(1,3,2)
imagesc(y(:,Late==1)',[0 7])
box off
subplot(1,3,3)
imagesc(z(:,Late==1)',[0 10])
box off

%%
clear all           
cd('T:\jan\Collabo Data\PFCpaperPreProcessed')
load('LickEvokedIndx.mat')
load AATC_PFC_Sua_Psth_1ms


AATC_Sua_Psth=AATC_PFC_Sua_Psth;

Bin=100;
AATC_Sua_Psth=squeeze(mean(reshape(AATC_Sua_Psth,Bin,size(AATC_Sua_Psth,1)/Bin,size(AATC_Sua_Psth,2),size(AATC_Sua_Psth,3))));
bc_Psths=AATC_Sua_Psth-mean(AATC_Sua_Psth(1:1000/Bin,:,:));
norm_bc_Psths=bc_Psths./max(bc_Psths);


% for i=1:size(bc_Psths,2)
%    bc_Psths(:,i,:) = smoothdata(bc_Psths(:,i,:),'gaussian',25);
% end


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
time=-1+Bin/1000:Bin/1000:4;

[Tr] = pca(PostExample(:,:,1)','NumComponents',3);
[Tur] = pca(PostExample(:,:,2)','NumComponents',3);
PFCexamplePopVectorCSplus=Tr;
PFCexamplePopVectorCSminus=Tur;

for i=1:3
   Tr(:,i) = smoothdata(Tr(:,i),'gaussian',7);
   Tur(:,i) = smoothdata(Tur(:,i),'gaussian',7);

end

PFCexamplePopVectorCSplusSmooth=Tr;
PFCexamplePopVectorCSminusSmooth=Tur;
cd('C:\Users\Jan\Desktop')

save('PFCexamplePopVectors','PFCexamplePopVectorCSplus','PFCexamplePopVectorCSminus','PFCexamplePopVectorCSplusSmooth','PFCexamplePopVectorCSminusSmooth')

plot3(Tr(base,1),Tr(base,2),Tr(base,3),'LineWidth',1,'Color',[0 0 1])
hold on
plot3(Tr(Tone,1),Tr(Tone,2),Tr(Tone,3),'LineWidth',1,'Color',[0 0 1])
hold on
plot3(Tr(Trace,1),Tr(Trace,2),Tr(Trace,3),'LineWidth',1,'Color',[0 0 1])
hold on
plot3(Tr(Reward,1),Tr(Reward,2),Tr(Reward,3),'LineWidth',1,'Color',[0 0 1])
hold on

plot3(Tur(base,1),Tur(base,2),Tur(base,3),'LineWidth',1,'Color',[1 0 0])
hold on
plot3(Tur(Tone,1),Tur(Tone,2),Tur(Tone,3),'LineWidth',1,'Color',[1 0 0])
hold on
plot3(Tur(Trace,1),Tur(Trace,2),Tur(Trace,3),'LineWidth',1,'Color',[1 0 0])
hold on
plot3(Tur(Reward,1),Tur(Reward,2),Tur(Reward,3),'LineWidth',1,'Color',[1 0 0])
hold on
xlabel('PCA1')
ylabel('PCA2')
zlabel('PCA3')
grid on
axis tight
hold on

figure()
hold on
    DistanceBC=(Distance-nanmean(Distance(1:(1000/Bin)-1,:)))./nanstd(Distance(1:(1000/Bin)-1,:));
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

[p,h]=signrank(nanmean(DistanceBC(30:40,Late==1)))


load SVM

figure()
n1{1}=1-allMcRb(LearnedCounter==1)
n1{2}=1-allMcRT(LearnedCounter==1)
names={'Baseline','Trace'};
% ttl={'SVM Classifier'};
ttl={''};
colors=[0 1 0;0 0 0;];
ylabel={'Percent Correct'}
style=([1 1 1;1 1 1]);
BeehivePlotRS(n1,names,ylabel,ttl,colors,style)
hold on
plot([.5 2.5],[0.5 0.5],':r','LineWidth',2)
set(gca,'FontSize',25);
set(gca,'LineWidth',5);


mean(n1{2})
[p,h]=signrank(allMcRT(LearnedCounter==1))

clearvars -except DistanceBC Late incl AvgLickCounter

cd('T:\jan\Collabo Data\PFCpaperPreProcessed')
    load Speed
    [ConFac]=ballCali()
  
  Bin=100*30;

  SpeedBinned=squeeze(nanmean(reshape(Speed,Bin,size(Speed,1)/Bin,2,size(Speed,3))))*(30000/Bin); 
  Speed=SpeedBinned(91:140,:,:);
  SpeedDiff=squeeze(Speed(:,1,:)-Speed(:,2,:));
  SpeedDiffZ=(SpeedDiff-mean(SpeedDiff(1:10,:)))./std(SpeedDiff(1:10,:));

LickPsth=AvgLickCounter(incl==1,:,1);
UnLickPsth=AvgLickCounter(incl==1,:,2);
Bin=100;
LickPsth=permute(LickPsth(:,:),[2,1,3]);
LickBinned=squeeze(sum(reshape(LickPsth,Bin,size(LickPsth,1)/Bin,size(LickPsth,2)))); 
UnLickPsth=permute(UnLickPsth(:,:),[2,1,3]);
UnLickBinned=squeeze(sum(reshape(UnLickPsth,Bin,size(UnLickPsth,1)/Bin,size(UnLickPsth,2)))); 
  
LickPsthInc=LickBinned(:,:);
UnLickPsthInc=UnLickBinned(:,:);
x=(LickPsthInc-UnLickPsthInc)
x=(x-nanmean(x(1:10,:)))./nanstd(x(1:10,:));
y=(DistanceBC(:,:)-(mean(DistanceBC(1:10,:))))./std(DistanceBC(1:10,:));
z=SpeedDiffZ(:,incl==1);
% 
% for i=1:length(find(incl==1))
% [h,p]=corr(x(10:40,i),y(10:40,i))
% Pcounter(i)=p;
% end
% SigCorPopLick=length(find(Pcounter<.05))

x(x==inf|x==-inf)=1;
x(isnan(x)==1)=1;
[h,p]=corrcoef(nanmean(x(10:40,Late==1)),nanmean(y(10:40,Late==1)))


% for i=1:length(find(incl==1))
% [h,p]=corr(x(10:40,i),z(10:40,i))
% Pcounter(i)=p;
% end
% SigCorPopSpeed=length(find(Pcounter<.05))

z(z==inf)=0
z(isnan(z)==1)=0
[h,p]=corrcoef(nanmean(z(10:40,Late==1)),nanmean(y(10:40,Late==1)))

[v,idx]=sort(Late);
subplot(1,3,1)
imagesc(x(:,Late==1)',[0 25])
box off
subplot(1,3,2)
imagesc(y(:,Late==1)',[0 7])
box off
subplot(1,3,3)
imagesc(z(:,Late==1)',[0 10])
box off


plot3(HPCexamplePopVectorCSminus(:,1),HPCexamplePopVectorCSminus(:,2),HPCexamplePopVectorCSminus(:,3),'LineWidth',1,'Color',[0 0 1])
