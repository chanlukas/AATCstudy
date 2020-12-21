clear all

%% Load and Combine All Event Data from HPConly Experiment
cd('T:\jan\Collabo Data')
load AATCR2hpcpfcAllBehavior.mat
Events={AATCR2Behavior.Events};
Animals={AATCR2Behavior.Animal};
Sessions={AATCR2Behavior.Session};

[ResponseRates]=HPCpaper_Pre_BehaviorLoop2R(Events,Animals,Sessions);



%% Plot Average Lick Responses


for i=1:5
ResponseRates1=ResponseRates{i};
ResponseRatesTraces(i,:,:)=squeeze(mean(ResponseRates1(1:10,:,30:40),3))-squeeze(mean(ResponseRates1(1:10,:,1:9),3));

end

ResponseRatesTraces([1,3],:,[1,2,3])=ResponseRatesTraces([1,3],:,[3,2,1])
figure()
pp(1,:)=nanmean(ResponseRatesTraces(:,:,1));
errBar=repmat(std(ResponseRatesTraces(:,:,1))/sqrt(size(ResponseRatesTraces,1)),2,1);
shadePlot2(1:10,pp(1,:),errBar,[0 0 1],.3)
hold on
box off
axis tight

pp(2,:)=nanmean(ResponseRatesTraces(:,:,2));
errBar=repmat(std(ResponseRatesTraces(:,:,2))/sqrt(size(ResponseRatesTraces,1)),2,1);
shadePlot2(1:10,pp(2,:),errBar,[1 0 0],.3)
hold on
box off
axis tight

pp(3,:)=nanmean(ResponseRatesTraces(:,:,3));
errBar=repmat(std(ResponseRatesTraces(:,:,3))/sqrt(size(ResponseRatesTraces,1)),2,1);
shadePlot2(1:10,pp(3,:),errBar,[0 0 1],.3)
hold on
box off
axis tight
ylabel('{\Delta} LickRate Trace Period [hz]')
xlabel('Training Days')

