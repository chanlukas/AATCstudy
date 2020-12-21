

function [PsthSua,EventMatrix,EventType]=HPCpaper_SuaPsth(SpikeTS,ClusterID,GoodClusters,EventTS,EventType,RecordingLength,Pre,Post,Bin)

for ii=1:length(GoodClusters)
   
    TScell(1:RecordingLength)=zeros(1,1,'int8');
    TScell(SpikeTS(find(ClusterID==GoodClusters(ii))))=ones(1,1,'int8');
        c=1;
        for i=1:length(EventTS)
           
        if EventTS(i)-((Pre*30000)-1)>0&&EventTS(i)+Post*30000<size(TScell,2)      %filtering out events for which the pre post cut off would violate recording bounds    

        EventMatrix(c,:,ii)=int8(TScell((EventTS(i)-((Pre*30000)-1)):(EventTS(i)+Post*(30000)))');   
        Eventgood(i)=1;
        c=c+1;
        else
        Eventgood(i)=0;
        end
        
        end
clearvars TScell
end
 
idx=unique(EventType);
if size(Eventgood,1)~=size(EventType,1)
Eventgood=Eventgood';
end  
EventType=EventType(Eventgood==1);
EventTS=EventTS(Eventgood==1);
for i = 1: length(idx)
EventMatrixMeans=squeeze(mean(EventMatrix(find(EventType==idx(i)),:,:),1)); 
if size(EventMatrixMeans,1)==1
   EventMatrixMeans=EventMatrixMeans';
end
PsthSua(:,:,i)=squeeze(sum(reshape(EventMatrixMeans,Bin*30,size(EventMatrixMeans,1)/(Bin*30),length(GoodClusters)))).*(1000/Bin);
StdErSua(:,:,i)=squeeze(std(reshape(EventMatrixMeans,Bin*30,size(EventMatrixMeans,1)/(Bin*30),length(GoodClusters))))/sqrt(size(EventMatrixMeans,2));
clearvars EventMatrixMeans
end
PsthSua=squeeze(PsthSua);

% % if Figure==1       
% for i=1:25%size(PsthSua,2);
%     if exist('Shank','var')==1
%     if Shank(i) ==1
%      AxisColor=[0 1 1];
%     elseif Shank(i) ==2
%      AxisColor=[1 1 0];
%     end
%     end
%  figure() 
% 
% set(0,'DefaultFigureWindowStyle','docked')
%     for ii=1:length(idx)
%        
% TS1=EventTS(EventType==idx(ii));
% 
% TS2(1:RecordingLength)=0;
% TS2=SpikeTS(find(ClusterID==GoodClusters(i)));
% 
% ts1=double(reshape(TS1,1,length(TS1)));
% ts2=double(reshape(TS2,1,length(TS2)));
% lags=[-Pre*30000,Post*30000];
% timestamps=[];
% trial=[];
% 
% for n = 1:length(ts1)
%     index=find(ts2>[ts1(n)+lags(1)]&ts2<[ts1(n)+lags(2)]);
%     timestamps=[timestamps,ts2(index)+(-ts1(n))];
%     trial=[trial,ones(1,length(index))*n];
% end
% if ii==1
%     color=[0 0 1];
% %     title1=('Rewarded Sound');
% elseif ii==2
%     color=[1 0 0];
% %     title1=('Unrewarded Sound');
% elseif ii==3
%     color=[0 0 .75];
% %     title1=('Rewarded Sound 2');
% end
% subplot(length(idx)+1,1,ii)
% % title(title1)
% hold on
% h=plot(timestamps/30000,trial,'s','markersize',2,'markerfacecolor',color,'markeredgecolor',color);
% axis([-Pre Post 0 length(EventTS(EventType==idx(ii)))])
% hold on
% plot([0 0],[0 max(trial)],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
% plot([2 2],[0 max(trial)],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
% plot([3 3],[0 max(trial)],'LineWidth',2,'LineStyle',':','Color',[0 0 0])
% ylabel('Trial')
% set(gca,'FontSize',25);
% set(gca,'LineWidth',5);
% 
% box off
% time=(-Pre+Bin/1000):Bin/1000:Post;
% subplot(length(idx)+1,1,length(idx)+1)
% errBar=repmat(StdErSua(:,i,ii),1,2)';
% pp=PsthSua(:,i,ii)';
% shadePlot(time,pp,errBar,color,.5)
% hold on
% plot([0 0],[max(max(pp))+max(max(errBar)) min(min(pp))-min(min(errBar))],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
% plot([2 2],[max(max(pp))+max(max(errBar)) min(min(pp))-min(min(errBar))],'LineWidth',2,'LineStyle','--','Color',[0 0 0])
% plot([3 3],[max(max(pp))+max(max(errBar)) min(min(pp))-min(min(errBar))],'LineWidth',2,'LineStyle',':','Color',[0 0 0])
% axis tight
% box off 
% ylabel('Firing Rate [Hz]')
% xlabel('Time [s]')
% set(gca,'FontSize',25);
% set(gca,'LineWidth',5);
% set(gca,'color')
%     end
% end
% end
end

