 
% Evaluates Single Cell Responses to AATC Task
clear all
%% SessionFinder
[Animals1,Sessions1]=AATC_SessionFinder()
[Animals2,Sessions2]=AATC_PFC_SessionFinder();
[Animals3,Sessions3]=AATC_HPC_PFC_SessionFinder();
[Animals4,Sessions4]=AATCR2_SessionFinder()
Animals=cat(2,Animals1,Animals2,Animals3,Animals4);
Sessions=cat(2,Sessions1,Sessions2,Sessions3,Sessions4);
Exp=cat(2,repmat(1,1,length(Animals1)),repmat(2,1,length(Animals2)),repmat(3,1,length(Animals3)),repmat(4,1,length(Animals4)));
%% Initialize Vars    
AATC_Sua_Psth=[];
wfs=[];
includeWF=[];
SpikeV2P=[];
spikeWidth=[];
Depth=[];
FRperiodsCounter=[];
SessionCounter=1;
CellCounter=0;

%% Main Loop
Speed=[];
c=1;
for a=find(ismember(Exp,[1,2,3]));
   
    for s= 1:size(Sessions{a},2)
    a
    s
    if Exp(a)==1
    cd(['T:\jan\HeadFixed Data\',Animals{a},'\',Sessions{a}{s},'\spikes'])
    elseif Exp(a)==2
    cd(['T:\jan\HeadFixed Data\HPC-PFC\',Animals{a},'\',Sessions{a}{s}])
    elseif Exp(a)==3
    cd(['T:\jan\HeadFixed Data\HPC128-PFC128\',Animals{a},'\',Sessions{a}{s}])
    elseif Exp(a)==4
    cd(['T:\jan\HeadFixed Data\AATCR2_HPC128-PFC128\',Animals{a},'\',Sessions{a}{s}])
    end
        
            try
        [SpeedPsth]= AATC_SpeedPsth([1,2],10,10);
        catch
        SpeedPsth(1:600000,1:2)=NaN;
        end

        Speed=cat(3,Speed,SpeedPsth);
   
        filename=dir('*rec*.txt');
Bpre=1999;
Bpost=4000;
Bbin=100;
BFig=0;
try
[RewardMatrix,unRewardMatrix]=AATCex(filename.name,Bpre,Bpost,Bbin,BFig); 

LickRate=mean(sum(RewardMatrix(:,4000:5000),2));
[Learned,LearnedP]=ttest2(mean(RewardMatrix(:,Bpre+2000:Bpre+2975),2),mean(unRewardMatrix(:,Bpre+2000:Bpre+2975),2));   
LearnedCounter(c)=Learned;
LearnedCounterP(c)=LearnedP;
catch
LearnedCounter(c)=NaN;
LearnedCounterP(c)=NaN;
end
c=c+1;
end
end


 
  cd('T:\jan\Collabo Data\PFCpaperPreProcessed')

   save AllSpeedHPCPFC Speed LearnedCounter LearnedCounterP
  
  [ConFac]=ballCali()
  
  Bin=100*30;

  SpeedBinned=squeeze(nanmean(reshape(Speed,Bin,size(Speed,1)/Bin,2,size(Speed,3))))*(30000/Bin); 
  
                idx1=find(LearnedCounter==1)
                idx2=find(LearnedCounter==1)
                figure()
                time=-9990:100:10000;
                PreR=squeeze(nanmean(SpeedBinned(:,1,idx1),3)*ConFac) ;              
                erPreR=repmat(nanstd(SpeedBinned(:,1,idx1),[],3)*ConFac,1,2);
                shadePlot(time,PreR',erPreR',[0 0 1],.5)
                plot([0 0],[-.5 2],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
                plot([2000 2000],[-.5 2],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
                plot([3000 3000],[-.5 2],'LineWidth',2,'LineStyle',':','Color',[0 0 0])
                  PreR=squeeze(nanmean(SpeedBinned(:,2,idx2),3)*ConFac) ;              
                erPreR=repmat(nanstd(SpeedBinned(:,2,idx2),[],3)*ConFac,1,2);
                shadePlot(time,PreR',erPreR',[1 0 0],.5)
                box off
                axis tight
                xlabel('Time [ms]')
                ylabel('Speed [cm/s]')      
                set(gca,'FontSize',25);
                set(gca,'LineWidth',5);
                