clear all

%% Load Data from HPConly Experiment
cd('T:\jan\Collabo Data')
load AATCstrucHPConly.mat
SuaGLM=[];
c=1;
for i=1:length(AATC_HPC)
    i
   SpikeTS=AATC_HPC(i).Spikes{1,1};
   ClusterID=AATC_HPC(i).Spikes{2,1};
   GoodClusters=unique(ClusterID);
   
   EventTypeAll=AATC_HPC(i).Events{2,1};
   EventTSAll=AATC_HPC(i).Events{1,1};
   
   EventType=EventTypeAll(ismember(EventTypeAll,[1,2]));
   EventTS=EventTSAll(ismember(EventTypeAll,[1,2]));
   
   RecordingLength=AATC_HPC(i).RecordingLength;
   
   Pre=1;
   Post=4;
   Bin=1;
   
   [PsthSua,EventMatrix,EventType]=HPCpaper_Pre_SuaPsth(SpikeTS,ClusterID,GoodClusters,EventTS,EventType,RecordingLength,Pre,Post,Bin);
    
   LickTS=EventTSAll(ismember(EventTypeAll,[3]));
   CSTS=EventTSAll(ismember(EventTypeAll,[1]));

   TimeFromLastLick=5*30000;
TimeFromLastCS=3*30000;
NextSecs=1*30000;
for ii=1:length(LickTS)

    if length(intersect(LickTS(ii)-TimeFromLastLick:LickTS(ii)-1,LickTS))==0&length(intersect(LickTS(ii)-TimeFromLastCS:LickTS(ii)-1,CSTS))==1 ...
       & length(intersect(LickTS(ii):LickTS(ii)+NextSecs,LickTS))>=3;
        LickBoutThres(ii)=1;
    else
        LickBoutThres(ii)=0;
    end
end

% figure()
LickTS=LickTS(LickBoutThres==1);

	[LickMatrixOG,EventIdxTimeOK]=HPCpaper_Pre_LickPsth(1*30000,4*30000,EventTS,EventType,LickTS,RecordingLength);
    bin=500*30;
    LickPredic=squeeze(mean(reshape(LickMatrixOG,size(LickMatrixOG,1),bin,size(LickMatrixOG,2)/bin),2))';

    
    
Rsound=zeros(size(find(EventIdxTimeOK==1),2),150000);
Rsound(EventType==1,30000:45000)=1;
Rsound=squeeze(mean(reshape(Rsound,size(LickMatrixOG,1),bin,size(LickMatrixOG,2)/bin),2))';

Base=zeros(size(find(EventIdxTimeOK==1),2),150000);
Base(:,1:30000)=1;
Base=squeeze(mean(reshape(Base,size(LickMatrixOG,1),bin,size(LickMatrixOG,2)/bin),2))';

URsound=zeros(size(find(EventIdxTimeOK==1),2),150000);
URsound(EventType==2,30000:45000)=1;
URsound=squeeze(mean(reshape(URsound,size(LickMatrixOG,1),bin,size(LickMatrixOG,2)/bin),2))';

Rsoundl=zeros(size(find(EventIdxTimeOK==1),2),150000);
Rsoundl(EventType==1,45000:90000)=1;
Rsoundl=squeeze(mean(reshape(Rsoundl,size(LickMatrixOG,1),bin,size(LickMatrixOG,2)/bin),2))';

URsoundl=zeros(size(find(EventIdxTimeOK==1),2),150000);
URsoundl(EventType==2,45000:90000)=1;
URsoundl=squeeze(mean(reshape(URsoundl,size(LickMatrixOG,1),bin,size(LickMatrixOG,2)/bin),2))';

Reward=zeros(size(find(EventIdxTimeOK==1),2),150000);
Reward(EventType==1,120000:135000)=1;
Reward=squeeze(mean(reshape(Reward,size(LickMatrixOG,1),bin,size(LickMatrixOG,2)/bin),2))';

Delay=zeros(size(find(EventIdxTimeOK==1),2),150000);
Delay(EventType==1,90000:120000)=1;
Delay=squeeze(mean(reshape(Delay,size(LickMatrixOG,1),bin,size(LickMatrixOG,2)/bin),2))';

DelayU=zeros(size(find(EventIdxTimeOK==1),2),150000);
DelayU(EventType==2,90000:120000)=1;
DelayU=squeeze(mean(reshape(DelayU,size(LickMatrixOG,1),bin,size(LickMatrixOG,2)/bin),2))';


cc=1;
for i=1:length(GoodClusters)
    i
cellFR=EventMatrix(:,:,i);
cellFR=squeeze(sum(reshape(cellFR,size(cellFR,1),bin,size(LickMatrixOG,2)/bin),2))';
glmData=table(Base(:),Rsound(:),Rsoundl(:),URsound(:),URsoundl(:),Reward(:),Delay(:),DelayU(:),LickPredic(:),cellFR(:));
mdl{cc} = fitglm(glmData,'Distribution','poisson');
cc=cc+1;

msg=lastwarn;
if size(msg,2)==0
BadScale(c)=0;
elseif size(msg,2)>0
BadScale(c)=1;    
end
    warning('')
c=c+1;
end

SuaGLM=[SuaGLM,mdl];
clearvars -except c i SuaGLM AATC_HPC BadScale
end

clearvars -except  c SuaGLM BadScale

load  AATCstrucHPCPFC.mat

for i=1:length(AATC_HPC_PFC)
    i
   SpikeTS=AATC_HPC_PFC(i).Spikes{1,1};
   ClusterID=AATC_HPC_PFC(i).Spikes{2,1};
   ShankS=cell2mat(AATC_HPC_PFC(i).Area);
   GoodClusters=unique(ClusterID);
   GoodClusters=GoodClusters(ShankS==1);
   
   EventTypeAll=AATC_HPC_PFC(i).Events{2,1};
   EventTSAll=AATC_HPC_PFC(i).Events{1,1};
   
   EventType=EventTypeAll(ismember(EventTypeAll,[1,2]));
   EventTS=EventTSAll(ismember(EventTypeAll,[1,2]));
   
   RecordingLength=AATC_HPC_PFC(i).RecordingLength;
   
   Pre=1;
   Post=4;
   Bin=1;
   
   [PsthSua,EventMatrix,EventType]=HPCpaper_Pre_SuaPsth(SpikeTS,ClusterID,GoodClusters,EventTS,EventType,RecordingLength,Pre,Post,Bin);

       LickTS=EventTSAll(ismember(EventTypeAll,[3]));
   CSTS=EventTSAll(ismember(EventTypeAll,[1]));

   TimeFromLastLick=5*30000;
TimeFromLastCS=3*30000;
NextSecs=1*30000;
for ii=1:length(LickTS)

    if length(intersect(LickTS(ii)-TimeFromLastLick:LickTS(ii)-1,LickTS))==0&length(intersect(LickTS(ii)-TimeFromLastCS:LickTS(ii)-1,CSTS))==1 ...
       & length(intersect(LickTS(ii):LickTS(ii)+NextSecs,LickTS))>=3;
        LickBoutThres(ii)=1;
    else
        LickBoutThres(ii)=0;
    end
end

% figure()
LickTS=LickTS(LickBoutThres==1);

	[LickMatrixOG,EventIdxTimeOK]=HPCpaper_Pre_LickPsth(1*30000,4*30000,EventTS,EventType,LickTS,RecordingLength);
    bin=500*30;
    LickPredic=squeeze(mean(reshape(LickMatrixOG,size(LickMatrixOG,1),bin,size(LickMatrixOG,2)/bin),2))';

    
    
Rsound=zeros(size(find(EventIdxTimeOK==1),2),150000);
Rsound(EventType==1,30000:45000)=1;
Rsound=squeeze(mean(reshape(Rsound,size(LickMatrixOG,1),bin,size(LickMatrixOG,2)/bin),2))';

Base=zeros(size(find(EventIdxTimeOK==1),2),150000);
Base(:,1:30000)=1;
Base=squeeze(mean(reshape(Base,size(LickMatrixOG,1),bin,size(LickMatrixOG,2)/bin),2))';

URsound=zeros(size(find(EventIdxTimeOK==1),2),150000);
URsound(EventType==2,30000:45000)=1;
URsound=squeeze(mean(reshape(URsound,size(LickMatrixOG,1),bin,size(LickMatrixOG,2)/bin),2))';

Rsoundl=zeros(size(find(EventIdxTimeOK==1),2),150000);
Rsoundl(EventType==1,45000:90000)=1;
Rsoundl=squeeze(mean(reshape(Rsoundl,size(LickMatrixOG,1),bin,size(LickMatrixOG,2)/bin),2))';

URsoundl=zeros(size(find(EventIdxTimeOK==1),2),150000);
URsoundl(EventType==2,45000:90000)=1;
URsoundl=squeeze(mean(reshape(URsoundl,size(LickMatrixOG,1),bin,size(LickMatrixOG,2)/bin),2))';

Reward=zeros(size(find(EventIdxTimeOK==1),2),150000);
Reward(EventType==1,120000:135000)=1;
Reward=squeeze(mean(reshape(Reward,size(LickMatrixOG,1),bin,size(LickMatrixOG,2)/bin),2))';

Delay=zeros(size(find(EventIdxTimeOK==1),2),150000);
Delay(EventType==1,90000:120000)=1;
Delay=squeeze(mean(reshape(Delay,size(LickMatrixOG,1),bin,size(LickMatrixOG,2)/bin),2))';

DelayU=zeros(size(find(EventIdxTimeOK==1),2),150000);
DelayU(EventType==2,90000:120000)=1;
DelayU=squeeze(mean(reshape(DelayU,size(LickMatrixOG,1),bin,size(LickMatrixOG,2)/bin),2))';


cc=1;
for i=1:length(GoodClusters)
    
cellFR=EventMatrix(:,:,i);
cellFR=squeeze(mean(reshape(cellFR,size(cellFR,1),bin,size(LickMatrixOG,2)/bin),2))';
glmData=table(Base(:),Rsound(:),Rsoundl(:),URsound(:),URsoundl(:),Reward(:),Delay(:),DelayU(:),LickPredic(:),cellFR(:));
mdl{cc} = fitglm(glmData,'Distribution','poisson');
cc=cc+1;

msg=lastwarn;
if size(msg,2)==0
BadScale(c)=0;
elseif size(msg,2)>0
BadScale(c)=1;    
end
    warning('')
c=c+1;
end

SuaGLM=[SuaGLM,mdl];
clearvars -except c i SuaGLM AATC_HPC_PFC BadScale

end


cd('T:\jan\Collabo Data\HPCpaperPreProcessed')
save('GLM.mat','SuaGLM','BadScale','-v7.3')

for i=1:size(SuaGLM,2)
SuaGLMpValues(i,:)=table2array(SuaGLM{1,i}.Coefficients(2:end,4));  
SuaGLMDirection(i,:)=table2array(SuaGLM{1,i}.Coefficients(2:end,1)); 
end    

save('SuaGLMpValues.mat','SuaGLMpValues','BadScale','SuaGLMDirection')
