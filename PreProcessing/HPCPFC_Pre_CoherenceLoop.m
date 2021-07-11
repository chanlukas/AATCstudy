  clear all

%% Load Data from HPConly Experiment
cd('/Users/Jan/Documents/PFC Paper/AATCdata')
load AATCstrucHPCPFC.mat
load AATChpcpfcAllBehavior.mat
AATC_Sua_Psth=[];
wfs=[];
includeWF=[];
SpikeV2P=[];
spikeWidth=[];
Depth=[];
FRperiodsCounter=[];
SessionCounter=1;
CellCounter=0;
CoherenceR=[];
CoherenceU=[];
%Exp1
for i=1:length(AATChpcpfcAllBehavior)  
BehaviorAnimals(i)=AATChpcpfcAllBehavior(i).Animal;
BehaviorSession(i)=AATChpcpfcAllBehavior(i).Session;
end
counter=1;
for i=1:length(AATC_HPC_PFC)
    
Learned(counter)=AATC_HPC_PFC(i).Learned(1);
LearnedP(counter)=AATC_HPC_PFC(i).Learned(2);
Session(counter)=cell2mat(AATC_HPC_PFC(i).Session(1));
Exp(counter)=cell2mat(AATC_HPC_PFC(i).Exp);
counter=counter+1;

Events=[1,2];
EventTS=AATC_HPC_PFC(i).Events{1,1}(ismember(AATC_HPC_PFC(i).Events{2,1},Events));
EventType=AATC_HPC_PFC(i).Events{2,1}(ismember(AATC_HPC_PFC(i).Events{2,1},Events));
unique(EventType)
pre=5*30000;
post=10*30000;
Eventgood=[];
 c=1;
 
CA1_Data=AATC_HPC_PFC(i).CA1_LFP;
PFC_Data=AATC_HPC_PFC(i).PFC_LFP;
    for ii =1 :size (EventTS,1)
        
        if EventTS(ii)-pre >0 && EventTS(ii)+post<size(CA1_Data,1) 
            startPoint=EventTS(ii)-pre;
            endPoint=EventTS(ii)+post;
            LFP_C(c,:)=downsample(double(CA1_Data(startPoint:endPoint)),100);
            LFP_P(c,:)=downsample(double(PFC_Data(startPoint:endPoint)),100);    
%             LFP_Cs(c,:)=smoothdata(downsample(double(CA1_Data(startPoint:endPoint)),100),'gaussian',20);
%             LFP_Ps(c,:)=smoothdata(downsample(double(PFC_Data(startPoint:endPoint)),100),'gaussian',20);

            c=c+1;
            Eventgood(ii)=1;
        else            
            Eventgood(ii)=0;
        end
    end
    
% NaNMatrix=zeros(size(LFP_C));
% for iii=1:size(LFP_C,1)
% meanHF=abs(diff(LFP_C(iii,:)));
% Threshold=500;
% for ii=101:size(meanHF,2)-100
% if (max(meanHF(ii-100:ii))>Threshold | max(meanHF(ii:ii+50))>Threshold)
% NaNMatrix(iii,ii,:) =1;
% end
% end
% end
% 
% LFP_C(NaNMatrix==1)=randi([-150 150],1);%(length(find(NaNMatrix==1))));
% LFP_P(NaNMatrix==1)=randi([-150 150],1);%(length(find(NaNMatrix==1))));

% Cmatrix_C= LFP_C-LFP_Cs;
% Cmatrix_P= LFP_P-LFP_Ps;
% 
% LFP_C(NaNMatrix==1)= Cmatrix_C(NaNMatrix==1);
% LFP_P(NaNMatrix==1)= Cmatrix_P(NaNMatrix==1);

params.Fs=300;
params.tapers=[];
params.fpass=[1 150];
params.trialave=0;
movingwin=[2 .1];
EventType=EventType(Eventgood==1);
[CohR,phi,S12,S1,S2,t,f]=cohgramc(LFP_C(EventType==1,:)',LFP_P(EventType==1,:)',movingwin,params); 
[CohU,phi,S12,S1,S2,t,f]=cohgramc(LFP_C(EventType==2,:)',LFP_P(EventType==2,:)',movingwin,params); 

% x=squeeze(mean(CohR,2));
% idx=find(max(x(1:51,:))<.6)
% 
% x1=squeeze(mean(CohU,2));
% idx1=find(max(x1(1:51,:))<.6)

CoherenceR=cat(3,CoherenceR,squeeze(mean(CohR(:,:,:),3)));
CoherenceU=cat(3,CoherenceU,squeeze(mean(CohU(:,:,:),3)));

% CoherenceLickTrig=cat(3,CoherenceLickTrig,CoherenceLT);

clearvars Events EventTS EventType CohR CohU LFP_C LFP_P CoherenceL LFP_Cs LFP_Ps Cmatrix_C Cmatrix_P NaNmatrix
i
end

cd('/Users/Jan/Documents/PFC Paper/AATCdata/PreprocessedDATA_PFCpaper')
save('HPCPFC_coherence','CoherenceR','CoherenceU','t','f','-v7.3')


 

imagesc(t,f,squeeze(nanmean(CohR(:,:,:),3)))