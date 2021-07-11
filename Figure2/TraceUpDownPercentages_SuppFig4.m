
%% Analysis sketches
clear all      
cd('/Users/Jan/Documents/PFC Paper/AATCdata/PreprocessedDATA_PFCpaper')
load AATC_Sua_Psth_1ms
% load SuaGlMpValues.mat

idx=find(diff(CellPerSesCounter)~=1)
mean(CellPerSesCounter(idx))

load('LickEvokedIndx.mat')
%% Trace Down
Bin=25;
window=3000/Bin:4000/Bin;
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


Tresh=-1;
Basenanstd=squeeze(nanstd(bc_Psths(baseline,:,:)));
EvR=squeeze(mean(bc_Psths(window,:,:)));
EvokedPeaks=squeeze(mean(bc_Psths(window,:,:)));
Evokedup(1:size(AATC_Sua_Psth,2))=0;
Evokedup(find(EvR(:,1)<Basenanstd(:,1)*Tresh|EvR(:,2)<Basenanstd(:,2)*Tresh))=1;

TraceDownHPC=Evokedup;


Condition2=find(LearnedCounter==1&Evokedup==1);
Condition1=find(LearnedCounter==0&TrgDayCounter<3&Evokedup==1);

data{1}=EvokedPeaks(Condition1,1)';
data{2}=EvokedPeaks(Condition1,2)';
data{3}=EvokedPeaks(Condition2,1)';
data{4}=EvokedPeaks(Condition2,2)';

Reward=zeros(size([data{1},data{2},data{3},data{4}],2),1);
Reward([1:size(data{1},2),(size(data{1},2)*2)+1:size(data{1},2)*2+1+size(data{3},2)])=1;
Learned=zeros(size([data{1},data{2},data{3},data{4}],2),1);
Learned(1:size(data{1},2)+size(data{2},2))=1;
[h,p]=anovan([data{1},data{2},data{3},data{4}],{Reward;Learned},'model','interaction')

[h,p]=ttest2(data{2},data{4})


AllPostCells=length(find(LearnedCounter==1))
PostTraceDown=length(find(LearnedCounter==1&Evokedup==1))
PerPostTracedown=PostTraceDown/AllPostCells

%% Trace Up

BaselineSTD=nanstd(AATC_Sua_PsthBined(baseline,:,:));
BaselineSTD(find(BaselineSTD==0))=NaN;
BaselineMean=nanmean(AATC_Sua_PsthBined(baseline,:,:));


bc_Psths=(AATC_Sua_PsthBined-BaselineMean)./BaselineSTD; %in herz
for i=1:size(bc_Psths,2)
   bc_Psths(:,i,:) = smoothdata(bc_Psths(:,i,:),'gaussian',20);
end



Tresh=1;
Basenanstd=squeeze(nanstd(bc_Psths(baseline,:,:)));
EvR=squeeze(mean((bc_Psths(window,:,:))));
EvokedPeaks=squeeze(mean(bc_Psths(window,:,:)));
Evokedup(1:size(AATC_Sua_Psth,2))=0;
Evokedup(find(EvR(:,1)>Basenanstd(:,1)*Tresh|EvR(:,2)>Basenanstd(:,2)*Tresh))=1;
MeanBaselineX=squeeze(mean(AATC_Sua_PsthBined(baseline,:,1)));

TraceUpHPC=Evokedup;


Condition2=find(LearnedCounter==1&Evokedup==1);
Condition1=find(LearnedCounter==0&TrgDayCounter<3&Evokedup==1);

data{5}=EvokedPeaks(Condition1,1)';
data{6}=EvokedPeaks(Condition1,2)';
data{7}=EvokedPeaks(Condition2,1)';
data{8}=EvokedPeaks(Condition2,2)';

Reward=zeros(size([data{5},data{6},data{7},data{8}],2),1);
Reward([1:size(data{5},2),(size(data{5},2)*2)+1:size(data{5},2)*2+1+size(data{7},2)])=1;
Learned=zeros(size([data{5},data{6},data{7},data{8}],2),1);
Learned(1:size(data{5},2)+size(data{6},2))=1;
[h,p]=anovan([data{5},data{6},data{7},data{8}],{Reward;Learned},'model','interaction')

ttest2(data{2},data{4})


AllPostCells=length(find(LearnedCounter==1))
PostTraceUp=length(find(LearnedCounter==1&Evokedup==1))
PerPostTraceUp=PostTraceUp/AllPostCells

%% Analysis sketches
clearvars -except TraceUpHPC TraceDownHPC data         
load AATC_PFC_Sua_Psth_1ms

load('LickEvokedIndx.mat')
AATC_Sua_Psth=AATC_PFC_Sua_Psth;



%% Trace Down
Bin=25;
window=3000/Bin:4000/Bin;
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



Tresh=-1;
Basenanstd=squeeze(nanstd(bc_Psths(baseline,:,:)));
EvR=squeeze(mean(bc_Psths(window,:,:)));
EvokedPeaks=squeeze(mean(bc_Psths(window,:,:)));
Evokedup(1:size(AATC_Sua_Psth,2))=0;
Evokedup(find(EvR(:,1)<Basenanstd(:,1)*Tresh|EvR(:,2)<Basenanstd(:,2)*Tresh))=1;

TraceDownPFC=Evokedup;


Condition1=find(LearnedCounter==1&Evokedup==1&LickDownPFC==0&LickUpPFC==0);
Condition2=find(LearnedCounter==0&TrgDayCounter<3&Evokedup==1);

data{9}=EvokedPeaks(Condition1,1)';
data{10}=EvokedPeaks(Condition1,2)';
data{11}=EvokedPeaks(Condition2,1)';
data{12}=EvokedPeaks(Condition2,2)';

Reward=zeros(size([data{9},data{10},data{11},data{12}],2),1);
Reward([1:size(data{9},2),(size(data{9},2)*2)+1:size(data{9},2)*2+1+size(data{11},2)])=1;
Learned=zeros(size([data{9},data{10},data{11},data{12}],2),1);
Learned(1:size(data{9},2)+size(data{10},2))=1;
[h,p]=anovan([data{9},data{10},data{11},data{12}],{Reward;Learned},'model','interaction')


AllPostCells=length(find(LearnedCounter==1))
PostTraceDown=length(find(LearnedCounter==1&Evokedup==1&LickDownPFC==0&LickUpPFC==0))
PerPostTracedown=PostTraceDown/AllPostCells
%% Trace Up
%Subtract Baseline and Normalizse

BaselineSTD=nanstd(AATC_Sua_PsthBined(baseline,:,:));
BaselineSTD(find(BaselineSTD==0))=NaN;
BaselineMean=nanmean(AATC_Sua_PsthBined(baseline,:,:));


bc_Psths=(AATC_Sua_PsthBined-BaselineMean)./BaselineSTD; %in herz
for i=1:size(bc_Psths,2)
   bc_Psths(:,i,:) = smoothdata(bc_Psths(:,i,:),'gaussian',20);
end



Tresh=1;
Basenanstd=squeeze(nanstd(bc_Psths(baseline,:,:)));
EvR=squeeze(mean((bc_Psths(window,:,:))));
EvokedPeaks=squeeze(mean(bc_Psths(window,:,:)));
Evokedup(1:size(AATC_Sua_Psth,2))=0;
Evokedup(find(EvR(:,1)>Basenanstd(:,1)*Tresh|EvR(:,2)>Basenanstd(:,2)*Tresh))=1;
MeanBaselineX=squeeze(mean(AATC_Sua_PsthBined(baseline,:,1)));

TraceUpPFC=Evokedup;


Condition2=find(LearnedCounter==1&Evokedup==1&LickDownPFC==0&LickUpPFC==0);
Condition1=find(LearnedCounter==0&TrgDayCounter<3&Evokedup==1);

data{13}=EvokedPeaks(Condition1,1)';
data{14}=EvokedPeaks(Condition1,2)';
data{15}=EvokedPeaks(Condition2,1)';
data{16}=EvokedPeaks(Condition2,2)';

Reward=zeros(size([data{13},data{14},data{15},data{16}],2),1);
Reward([1:size(data{13},2),(size(data{13},2)*2)+1:size(data{13},2)*2+1+size(data{14},2)])=1;
Learned=zeros(size([data{13},data{14},data{15},data{16}],2),1);
Learned(1:size(data{13},2)+size(data{14},2))=1;
[h,p]=anovan([data{13},data{14},data{15},data{16}],{Reward;Learned},'model','interaction')



[H,P]=ttest2(data{14},data{16})

AllPostCells=length(find(LearnedCounter==1))
PostTraceUp=length(find(LearnedCounter==1&Evokedup==1&LickDownPFC==0&LickUpPFC==0))
PerPostTraceUp=PostTraceUp/AllPostCells


y=[nanmean(data{1}),nanmean(data{2}),nanmean(data{3}),nanmean(data{4}),...
    nanmean(data{5}),nanmean(data{6}),nanmean(data{7}),nanmean(data{8}),...
    nanmean(data{9}),nanmean(data{10}),nanmean(data{11}),nanmean(data{12}),...
    nanmean(data{13}),nanmean(data{14}),nanmean(data{15}),nanmean(data{16})]

e=[nanstd(data{1})/sqrt(length(data{1})),...
nanstd(data{2})/sqrt(length(data{2}))....
nanstd(data{3})/sqrt(length(data{3})),...
nanstd(data{4})/sqrt(length(data{4})),...
nanstd(data{5})/sqrt(length(data{5})),...
nanstd(data{6})/sqrt(length(data{6}))....
nanstd(data{7})/sqrt(length(data{7})),...
nanstd(data{8})/sqrt(length(data{8})),...
nanstd(data{9})/sqrt(length(data{9})),...
nanstd(data{10})/sqrt(length(data{10}))....
nanstd(data{11})/sqrt(length(data{11})),...
nanstd(data{12})/sqrt(length(data{12})),...
nanstd(data{13})/sqrt(length(data{13})),...
nanstd(data{14})/sqrt(length(data{14}))....
nanstd(data{15})/sqrt(length(data{15})),...
nanstd(data{16})/sqrt(length(data{16}))];


p=[ ]
c=[0 0 1;1 0 0]
superbar(y,'E',e,'P',p,'BarFaceColor',c)
names={'CS+','CS-'}
set(gca,'Xtick',0:1:3)
set(gca,'xticklabel',names)


%% Analysis sketches
clear all      
cd('/Users/Jan/Documents/PFC Paper/AATCdata/PreprocessedDATA_PFCpaper')
load AATC_Sua_Psth_1ms
load('LickEvokedIndx.mat')
%% Trace Down
Bin=25;
window=3000/Bin:4000/Bin;
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


Tresh=-1;
Basenanstd=squeeze(nanstd(bc_Psths(baseline,:,:)));
EvR=squeeze(mean(bc_Psths(window,:,:)));

EvokedupR(1:size(AATC_Sua_Psth,2))=0;
EvokedupR(find(EvR(:,1)<Basenanstd(:,1)*Tresh))=1;

EvokedupU(1:size(AATC_Sua_Psth,2))=0;
EvokedupU(find(EvR(:,2)<Basenanstd(:,2)*Tresh))=1;


LengthPre=length(find(LearnedCounter==0&TrgDayCounter<3))
LengthPost=length(find(LearnedCounter==1))


data{1}=length(find(LearnedCounter==0&TrgDayCounter<3&EvokedupR==1))/LengthPre;
data{2}=length(find(LearnedCounter==0&TrgDayCounter<3&EvokedupU==1))/LengthPre;
data{3}=length(find(LearnedCounter==1&EvokedupR==1))/LengthPost;
data{4}=length(find(LearnedCounter==1&EvokedupU==1))/LengthPost;


%% Trace Up

BaselineSTD=nanstd(AATC_Sua_PsthBined(baseline,:,:));
BaselineSTD(find(BaselineSTD==0))=NaN;
BaselineMean=nanmean(AATC_Sua_PsthBined(baseline,:,:));


bc_Psths=(AATC_Sua_PsthBined-BaselineMean)./BaselineSTD; %in herz
for i=1:size(bc_Psths,2)
   bc_Psths(:,i,:) = smoothdata(bc_Psths(:,i,:),'gaussian',20);
end



Tresh=1;
Basenanstd=squeeze(nanstd(bc_Psths(baseline,:,:)));
EvR=squeeze(mean((bc_Psths(window,:,:))));

EvokedupR(1:size(AATC_Sua_Psth,2))=0;
EvokedupR(find(EvR(:,1)>Basenanstd(:,1)*Tresh))=1;

EvokedupU(1:size(AATC_Sua_Psth,2))=0;
EvokedupU(find(EvR(:,2)>Basenanstd(:,2)*Tresh))=1;


LengthPre=length(find(LearnedCounter==0&TrgDayCounter<3))
LengthPost=length(find(LearnedCounter==1))


data{5}=length(find(LearnedCounter==0&TrgDayCounter<3&EvokedupR==1))/LengthPre;
data{6}=length(find(LearnedCounter==0&TrgDayCounter<3&EvokedupU==1))/LengthPre;
data{7}=length(find(LearnedCounter==1&EvokedupR==1))/LengthPost;
data{8}=length(find(LearnedCounter==1&EvokedupU==1))/LengthPost;
%% Analysis sketches
clearvars -except data         
load AATC_PFC_Sua_Psth_1ms

load('LickEvokedIndx.mat')
AATC_Sua_Psth=AATC_PFC_Sua_Psth;



%% Trace Down
Bin=25;
window=3000/Bin:4000/Bin;
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


Tresh=-1;
Basenanstd=squeeze(nanstd(bc_Psths(baseline,:,:)));
EvR=squeeze(mean(bc_Psths(window,:,:)));

EvokedupR(1:size(AATC_Sua_Psth,2))=0;
EvokedupR(find(EvR(:,1)<Basenanstd(:,1)*Tresh))=1;

EvokedupU(1:size(AATC_Sua_Psth,2))=0;
EvokedupU(find(EvR(:,2)<Basenanstd(:,2)*Tresh))=1;


LengthPre=length(find(LearnedCounter==0&TrgDayCounter<3))
LengthPost=length(find(LearnedCounter==1))


data{9}=length(find(LearnedCounter==0&TrgDayCounter<3&EvokedupR==1))/LengthPre;
data{10}=length(find(LearnedCounter==0&TrgDayCounter<3&EvokedupU==1))/LengthPre;
data{11}=length(find(LearnedCounter==1&EvokedupR==1))/LengthPost;
data{12}=length(find(LearnedCounter==1&EvokedupU==1))/LengthPost;


%% Trace Up

BaselineSTD=nanstd(AATC_Sua_PsthBined(baseline,:,:));
BaselineSTD(find(BaselineSTD==0))=NaN;
BaselineMean=nanmean(AATC_Sua_PsthBined(baseline,:,:));


bc_Psths=(AATC_Sua_PsthBined-BaselineMean)./BaselineSTD; %in herz
for i=1:size(bc_Psths,2)
   bc_Psths(:,i,:) = smoothdata(bc_Psths(:,i,:),'gaussian',20);
end



Tresh=1;
Basenanstd=squeeze(nanstd(bc_Psths(baseline,:,:)));
EvR=squeeze(mean((bc_Psths(window,:,:))));

EvokedupR(1:size(AATC_Sua_Psth,2))=0;
EvokedupR(find(EvR(:,1)>Basenanstd(:,1)*Tresh))=1;

EvokedupU(1:size(AATC_Sua_Psth,2))=0;
EvokedupU(find(EvR(:,2)>Basenanstd(:,2)*Tresh))=1;


LengthPre=length(find(LearnedCounter==0&TrgDayCounter<3))
LengthPost=length(find(LearnedCounter==1))


data{13}=length(find(LearnedCounter==0&TrgDayCounter<3&EvokedupR==1))/LengthPre;
data{14}=length(find(LearnedCounter==0&TrgDayCounter<3&EvokedupU==1))/LengthPre;
data{15}=length(find(LearnedCounter==1&EvokedupR==1))/LengthPost;
data{16}=length(find(LearnedCounter==1&EvokedupU==1))/LengthPost;

y=[nanmean(data{1}),nanmean(data{2}),nanmean(data{3}),nanmean(data{4}),...
    nanmean(data{5}),nanmean(data{6}),nanmean(data{7}),nanmean(data{8}),...
    nanmean(data{9}),nanmean(data{10}),nanmean(data{11}),nanmean(data{12}),...
    nanmean(data{13}),nanmean(data{14}),nanmean(data{15}),nanmean(data{16})]

e=[nanstd(data{1})/sqrt(length(data{1})),...
nanstd(data{2})/sqrt(length(data{2}))....
nanstd(data{3})/sqrt(length(data{3})),...
nanstd(data{4})/sqrt(length(data{4})),...
nanstd(data{5})/sqrt(length(data{5})),...
nanstd(data{6})/sqrt(length(data{6}))....
nanstd(data{7})/sqrt(length(data{7})),...
nanstd(data{8})/sqrt(length(data{8})),...
nanstd(data{9})/sqrt(length(data{9})),...
nanstd(data{10})/sqrt(length(data{10}))....
nanstd(data{11})/sqrt(length(data{11})),...
nanstd(data{12})/sqrt(length(data{12})),...
nanstd(data{13})/sqrt(length(data{13})),...
nanstd(data{14})/sqrt(length(data{14}))....
nanstd(data{15})/sqrt(length(data{15})),...
nanstd(data{16})/sqrt(length(data{16}))];

figure()
p=[ ]
c=[0 0 1;1 0 0]
superbar(y,'E',e,'P',p,'BarFaceColor',c)
names={'CS+','CS-'}
set(gca,'Xtick',0:1:3)
set(gca,'xticklabel',names)


d
