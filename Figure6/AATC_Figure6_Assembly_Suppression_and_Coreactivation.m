% % Code for reproducing Figure 6 and associated Supp. Figures
% % Requires Cell-Assembly-Detection package on path (see https://github.com/tortlab/Cell-Assembly-Detection)
% % Uses sub.m auxiliary function to plot.

clear
load('preprocessed_data', 'data', 'hasPF', 'hasCA1', 'learned', 'ses', 'exp', 'trainingday', 'animal',...
    'varnames', 'srate', 'datafiles')

%%

day=[];
animal=[];
for i=1:length(data)
    day(i)= data(i).TrainingDay{1}(1);
    animal(i)= data(i).Animal{1}; 
end

%%
nlick_all={};
for ises = 1:length(data)
    disp(['Performing session: ' num2str(ises)])
    animal(ises) = data(ises).Animal{1};
    spktimes = double(data(ises).Spikes{1})/srate;
    spkid = data(ises).Spikes{2};
    area = data(ises).Area{1};
    etime = data(ises).Events{1}/srate;
    etype = data(ises).Events{2};
    spkid = spkid-min(spkid)+1;
    classes = unique(spkid);

    %%
    idx1 = find(etype == 1);
    trialtype = 1*ones(length(idx1),1);
    idx2 = find(etype == 2);
    trialtype = [trialtype; 2*ones(length(idx2),1)];
    idx = [idx1; idx2];
    trialtime = etime(idx);
    trialtime(:,2) = etime(idx)+3;
    
    valid = 0;
%     
%     timebins = [-1 0;
%                 0 2;
%                 2 3;
%                 3 4];

    step=200/1000;
    timebins = -1:step:4;
    binwidth=step;
    preidx = find(timebins>-1 & timebins+binwidth<=0);
    traceidx = find(timebins>2 & timebins+binwidth<=3);
    
    nlick=[];
    nlick_isi=[];
    nlick_first=[];
    for itrial = 1:size(trialtime,1)
        valid = etime>=trialtime(itrial,1)-1 & etime<trialtime(itrial,2)+1;
        valid = valid & etype==3;
%         return
        etime_ = etime(valid);
        etime_ = etime_-trialtime(itrial,1);
        etype_ = etype(valid);
        
        for ibin=1:length(timebins)-1
            idx = find(etime_>timebins(ibin) & etime_<=timebins(ibin+1));
            nlick(itrial,ibin) = length(idx);
            
        end
        idx = find(etime_>0 & etime_<=2);
        nlick_isi(itrial,1) = mean(diff(etime_(idx)));
        
        nlick_first(itrial,1)=NaN;
        if ~isempty(idx)
            nlick_first(itrial,1) = min(etime_(idx));
        end
        
        idx = find(etime_>2 & etime_<=3);
        nlick_isi(itrial,2) = mean(diff(etime_(idx)));
        
        nlick_first(itrial,2)=NaN;
        if ~isempty(idx)
            nlick_first(itrial,2) = min(etime_(idx));
        end
        
        idx = find(etime_>3 & etime_<=4);
        nlick_isi(itrial,3) = mean(diff(etime_(idx)));
        nlick_first(itrial,3)=NaN;
        if ~isempty(idx)
            nlick_first(itrial,3) = min(etime_(idx));
        end
        
        
        %%
    end
    trialtype_all{ises} = trialtype;
    trialtime_all{ises} = trialtime;
    nlick_all{ises} = nlick;
    nlick_isi_all{ises} =nlick_isi;
    nlick_first_all{ises} =nlick_first;
end


idx_pre = find(timebins<-0.1);
idx_stim = find(timebins>=-0.1 & timebins<2);
idx_trace = find(timebins>=2 & timebins<3);
idx_rwd = find(timebins>=3);
%%

CA1=1;
PFC=2;

PRE = 1;
STIM = 2;
TRACE = 3;
RWD = 4;
ALL = 5;


timebins = -1:0.001:4;
periods={};
periods{PRE} = find(timebins<0);
periods{STIM} = find(timebins>=0 & timebins<0.5);
periods{TRACE} = find(timebins>=2 & timebins<3);
periods{RWD} = find(timebins>=3 & timebins<3.5);
periods{ALL} = find(timebins>=-Inf);


doplot = false;
valid_ses=[hasCA1', hasPF', hasCA1' & hasPF'];
% valid_ses={hasCA1, hasPF, hasCA1 & hasPF};
valid_neurons={CA1, PFC, []};
nsub = 20;
binwitdh = 20;

area_labels={'CA1', 'PFC', 'CA1 & PFC'};

%%

clear mean_assembly_all
mean_assembly_all={};

try
    load(['Assembly_patterns'])
    computeassemblies=0;
catch
    computeassemblies=1;
end

%%
try
    load(['Figure6_data'], 'area_labels',...
        'CA1','PFC','binwitdh','valid_neurons','valid_ses',...
        'mean_assembly_all','nneurons_all')
catch
    
    for ises = 1:size(valid_ses,1) 
        clear mean_assembly

        %seed = 654299090;
        %rng(seed)


        animal(ises) = data(ises).Animal{1};
        spktimes = double(data(ises).Spikes{1})/srate;
        spkid = data(ises).Spikes{2};
        area = data(ises).Area{1};
        etime = data(ises).Events{1}/srate;
        etype = data(ises).Events{2};
        spkid = spkid-min(spkid)+1;
        classes = unique(spkid);

        if isempty(data(ises).RipTS)
            continue
        end
        riptime = double(data(ises).RipTS{1})/srate;

        area_aux =[];
        area_aux(classes) =area;


        lessneuron=[];

        idx1 = find(etype == 1);
        trialtype = 1*ones(length(idx1),1);
        idx2 = find(etype == 2);
        trialtype = [trialtype; 2*ones(length(idx2),1)];
        idx = [idx1; idx2];
        trialtime = etime(idx);
        trialtime(:,2) = etime(idx)+3;

        timebins = -1:0.001:4;

        nneurons=[];
        for iarea=1:length(area_labels)-1
            nneurons(iarea) = sum(area_aux(unique(spkid))==valid_neurons{iarea});
        end

        if sum(nneurons==0)
            nneurons(iarea+1) = 0;
        else
            nneurons(iarea+1)= sum(nneurons(1:2));
        end

        %%
        for iarea = find(nneurons~=0)

            % No subsampling needed here
            nsubsample=1;
            for isample = 1:nsubsample
                area_=area;
                area_both=area;

                valid = 0;


                timebins = -1:0.001:4;
                periods={};
                periods{PRE} = find(timebins<0);
                periods{STIM} = find(timebins>=0 & timebins<1);
                periods{TRACE} = find(timebins>=2 & timebins<3);
                periods{RWD} = find(timebins>=3 & timebins<4);
                periods{ALL} = find(timebins>=-Inf);

                for isource=ALL
                %%
                timebins = -1:0.001:4;
                psth=zeros(length(classes),size(trialtime,1),length(timebins));
                psth_smooth=psth;
                for itrial = 1:size(trialtime,1)
                    valid = spktimes>=trialtime(itrial,1)-1 & spktimes<trialtime(itrial,2)+1;

                    if ~isempty(valid_neurons{iarea})
                        area_aux =[];
                        area_aux(classes) =area_;
                        valid = valid & (area_aux(spkid)==valid_neurons{iarea})';
                    else
                        area_aux =[];
                        area_aux(classes) =area_both;
                        valid = valid & (area_aux(spkid)>0)';
                    end
                    spktimes_ = spktimes(valid);
                    spktimes_ = spktimes_-trialtime(itrial,1);

                    spkid_ = spkid(valid);

                    for ineuron=1:length(classes)
                        idx = find(spkid_==classes(ineuron));
                        if isempty(idx)
                            continue
                        end
                        count = histc(spktimes_(idx),timebins);
                        psth(ineuron,itrial,:)=count;
                        psth_smooth(ineuron,itrial,:)=conv(count,rectwin(binwitdh),'same');
                    end
                end
                %%

                nlick = nlick_all{ises};
                lick_trace = mean(nlick(:,idx_trace),2);            
                ntrial = min(sum(trialtype==1),sum(trialtype==2));

                right1 = find(trialtype==1 & lick_trace > 0);
                wrong1 = find(trialtype==1 & lick_trace <= 0);
                wrong2 = find(trialtype==2 & lick_trace > 0);
                right2 = find(trialtype==2 & lick_trace <= 0);

                %%

                if computeassemblies==1
                    M = psth_smooth(:,:,periods{isource});
                    M = M(:,:,1:20:end);
                    M = reshape(M,size(M,1),[]);
                    P = assembly_patterns(M);
                else
                    P=patterns_all{iarea,ises,isource};
                end
                if isempty(P)
                    continue
                end

                psth_proj=zeros(size(P,2),size(psth_smooth,2),length(timebins));
                parfor itrial=1:size(psth_smooth,2)
                    psth_proj(:,itrial,:) = assembly_activity(P,squeeze(psth_smooth(:,itrial,:)));
                end            


                mean_assembly_stim = squeeze(nanmean(psth_proj,2));
                sem_assembly_stim = squeeze(nanstd(psth_proj,0,2)./...
                    sqrt(sum(~isnan(psth_proj),2)));

                mean_right_stim1 = squeeze(nanmean(psth_proj(:,right1,:),2));
                mean_wrong_stim1 = squeeze(nanmean(psth_proj(:,wrong1,:),2));

                mean_right_stim2 = squeeze(nanmean(psth_proj(:,right2,:),2));
                mean_wrong_stim2 = squeeze(nanmean(psth_proj(:,wrong2,:),2));


                %%
                timebins = -1:0.001:1;
                licktime = etime(etype==3);
                licktime(:,2) = licktime+1;
                clear psth psth_smooth
                psth = zeros(length(classes),size(licktime,1),length(timebins));
                psth_smooth = zeros(length(classes),size(licktime,1),length(timebins));
    %             tic
                for ilick = 1:size(licktime,1)
                    valid = spktimes>=licktime(ilick,1)-1 & spktimes<licktime(ilick,2);

                    if ~isempty(valid_neurons{iarea})
                        area_aux =[];
                        area_aux(classes) =area_;
                        valid = valid & (area_aux(spkid)==valid_neurons{iarea})';
                    else
                        area_aux =[];
                        area_aux(classes) =area_both;
                        valid = valid & (area_aux(spkid)>0)';
                    end
                    spktimes_ = spktimes(valid);
                    spktimes_ = spktimes_-licktime(ilick,1);

                    spkid_ = spkid(valid);

                    for ineuron=1:length(classes)
                        idx = find(spkid_==classes(ineuron));
                        if isempty(idx)
                            continue
                        end
                        count = histc(spktimes_(idx),timebins);
                        psth(ineuron,ilick,:)=count;
                        psth_smooth(ineuron,ilick,:)=conv(count,rectwin(binwitdh),'same');
                        if doplot
                            if ~isempty(idx)
                                plot(spktimes_(idx),ineuron,'.k')
                            end
                        end
                    end
                    %%
                end
    %             toc
                %%
                psth_proj_lick=zeros(size(P,2),size(psth_smooth,2),length(timebins));
                parfor ilick=1:size(psth_smooth,2)
                    psth_proj_lick(:,ilick,:) = assembly_activity(P,squeeze(psth_smooth(:,ilick,:)));
                end            


                mean_assembly_lick = squeeze(nanmean(psth_proj_lick,2));
                sem_assembly_lick = squeeze(nanstd(psth_proj_lick,0,2)./sqrt(sum(~isnan(psth_proj_lick),2)));



                %%


                timebins = -1:0.001:1;
                riptime = double(data(ises).RipTS{1})/srate;


                riptime=riptime(:);
                riptime(:,2) = riptime+1;
                rip_prevtrial = zeros(1,size(riptime,1));
                rip_nexttrial = zeros(1,size(riptime,1));
                clear psth psth_smooth
                psth = zeros(length(classes),size(riptime,1),length(timebins));
                psth_smooth = zeros(length(classes),size(riptime,1),length(timebins));

                for irip = 1:size(riptime,1)


                    valid = spktimes>=riptime(irip,1)-1 & spktimes<riptime(irip,2);


                    [v1,idx] = sort(trialtime(:,1));
                    [v2,iidx] = sort(idx);

                    aux = find(riptime(irip,1)<=trialtime(idx(2:end),1)-1 & riptime(irip,1)>trialtime(idx(1:end-1),2)+1);

                    if ~isempty(aux)
                        rip_prevtrial(irip) = find(iidx==aux);
                        rip_nexttrial(irip) = find(iidx==aux+1);
                    end


                    if ~isempty(valid_neurons{iarea})
                        area_aux =[];
                        area_aux(classes) =area_;
                        valid = valid & (area_aux(spkid)==valid_neurons{iarea})';
                    else
                        area_aux =[];
                        area_aux(classes) =area_both;
                        valid = valid & (area_aux(spkid)>0)';
                    end
                    spktimes_ = spktimes(valid);
                    spktimes_ = spktimes_-riptime(irip,1);

                    spkid_ = spkid(valid);

                    for ineuron=1:length(classes)
                        idx = find(spkid_==classes(ineuron));
                        if isempty(idx)
                            continue
                        end
                        count = histc(spktimes_(idx),timebins);
                        psth(ineuron,irip,:)=count;
                        psth_smooth(ineuron,irip,:)=conv(count,rectwin(binwitdh),'same');
                        if doplot
    %                         if ~isempty(idx)
    %                             plot(spktimes_(idx),ineuron,'.k')
    %                         end
                        end
                    end
                    %%
                end

                psth_proj_rip=zeros(size(P,2),size(psth_smooth,2),length(timebins));
                parfor irip=1:size(psth_smooth,2)
                    psth_proj_rip(:,irip,:) = assembly_activity(P,squeeze(psth_smooth(:,irip,:)));
                end                        

                mean_assembly_rip = squeeze(nanmean(psth_proj_rip,2));
                n_rip = [];
                idx = find(sum(right1==rip_prevtrial));
                mean_right1_rip = squeeze(nanmean(psth_proj_rip(:,idx,:),2));            
                n_rip = [n_rip length(idx)];

    %             idx = intersect(right2,rip_prevtrial);
                idx = find(sum(right2==rip_prevtrial));
                mean_right2_rip = squeeze(nanmean(psth_proj_rip(:,idx,:),2));
                n_rip = [n_rip length(idx)];


                idx = find(sum(wrong1==rip_prevtrial));
                mean_wrong1_rip = squeeze(nanmean(psth_proj_rip(:,idx,:),2));
                n_rip = [n_rip length(idx)];

                idx = find(sum(wrong2==rip_prevtrial));
                mean_wrong2_rip = squeeze(nanmean(psth_proj_rip(:,idx,:),2));
                n_rip = [n_rip length(idx)];


                sem_assembly_rip = squeeze(nanstd(psth_proj_rip,0,2)./sqrt(sum(~isnan(psth_proj_rip),2)));


                %%
                n_stim = [length(right1) length(right2) length(wrong1) length(wrong2)];

                mean_assembly{1,1} = mean_right_stim1;
                mean_assembly{2,1} = mean_right_stim2;
                mean_assembly{3,1} = mean_assembly_stim;
                mean_assembly{4,1} = mean_assembly_lick;
                mean_assembly{5,1} = mean_assembly_rip;

                mean_assembly{1,2} = mean_right_stim1;
                mean_assembly{2,2} = mean_right_stim2;
                mean_assembly{3,2} = mean_wrong_stim1;
                mean_assembly{4,2} = mean_wrong_stim2;
                mean_assembly{5,2} = n_stim;

                mean_assembly{1,3} = mean_right1_rip;
                mean_assembly{2,3} = mean_right2_rip;
                mean_assembly{3,3} = mean_wrong1_rip;
                mean_assembly{4,3} = mean_wrong2_rip;
                mean_assembly{5,3} = n_rip;

                mean_assembly_all{iarea,ises,isource} = mean_assembly;

    %             patterns_all{iarea,ises,isource} = P;

                end

            if isempty(lessneuron)
                nneurons_all(iarea,ises) = size(psth,1);
            else
                nneurons_all(iarea,ises) = lessneuron;
            end

            nass_all(iarea,ises) = size(P,2);

            ntrial = min(sum(trialtype==1),sum(trialtype==2));
            ntrials_all(iarea,ises) = ntrial;

            end     
        end

        %%
        save(['Figure6_data'], 'area_labels',...
            'CA1','PFC','binwitdh','valid_neurons','valid_ses',...
            'mean_assembly_all','nneurons_all','-v7.3')


    end
end

%% Plotting Figure 6
    
isource=ALL;
CSPLUS=1;
CSMINUS=2;
STIM=3;
LICK=4;
RIP = 5;

learn_labels = {'learning' 'overtrained'}
col = [.2 .7 .2;
    .2 .2 .7;
    .7 .2 .2];


col = [.9 .7 .7;
    .7 .7 .9;
    .7 .7 .7;
    .7 .2 .2;
    .2 .2 .7;
    0 0 0];

col2 = [.7 .7 .7;
    0 0 0];

figure(1)
clf

sub(4,3,4,1:2); cla
hold on
plot(0,nan,'color',col(1,:),'linewidth',2)
plot(0,nan,'color',col(2,:),'linewidth',2)
plot(0,nan,'color',col(3,:),'linewidth',2)
plot(0,nan,'color',col(4,:),'linewidth',2)
plot(0,nan,'color',col(5,:),'linewidth',2)
plot(0,nan,'color',col(6,:),'linewidth',2)

figure(2); clf

peak_area = [1 1 2];
maux_all={}
for iarea = 1:2
for iperiod = 1:2
    
    modulation_all=[];
    activation_all=[];
    mrip_all=[];

    for ilearn = 1:2
        
        %%
        idx = find(learned(:,1)==ilearn-1 & valid_ses(:,iarea));
        ises = idx;

        ises_area = find(valid_ses(:,iarea));

        m_all_rip_hi=[];
        m_all_rip_lo=[];
        m_all_rip_n=[];
        m_all_rip=[];
        m_all_stim1=[];
        m_all_stim2=[];
        m_all_lick=[];
        for ises_ = ises'

            mean_assembly = mean_assembly_all{iarea,ises_,5};

            if isempty(mean_assembly)
                continue
            end
            
            P = patterns_all{iarea,ises_,isource};
            
            ass_norm=[];
            for iass = 1:size(P,2)
                w=P(:,iass);
                T = w*w';
                T = T-diag(diag(T)); 
                ass_norm(iass) = ((w)'*T*(w));
            end
             
            mod = (mean_assembly{CSPLUS,1}-mean_assembly{CSMINUS,1});
            if size(mod,2)==1
                mod=mod';
            end
            
            
            timebins = -1:0.001:4;
            
            idx_stim = find(timebins>=0 & timebins<=2);
            idx_trace = find(timebins>=2 & timebins<=3);
            mod = bsxfun(@rdivide, mod, ass_norm');
            
            if iperiod==1
                [v, modhi] = max(mean(mod(:,idx_stim),2));
                [v, modlo] = min(mean(mod(:,idx_stim),2));
                [v, modn] = min(mean(abs(mod(:,idx_stim)),2));
            elseif iperiod==2
                [v, modhi] = max(mean(mod(:,idx_trace),2));
                [v, modlo] = min(mean(mod(:,idx_trace),2));
                [v, modn] = min(mean(abs(mod(:,idx_trace)),2));
            end
            
            %%
            
            m = mean_assembly{RIP,1};
            if size(m,2)==1
                m=m';
            end
            
            m = bsxfun(@rdivide, m, ass_norm');
            
            
            %% Get the assemblies with higher and less modulation.. 
            
            %% Correlate CS+ CS- diff with assembly reactivation!!
            m_all_rip_hi = cat(1,m(modhi,:),m_all_rip_hi);
            m_all_rip_lo = cat(1,m(modlo,:),m_all_rip_lo);
            m_all_rip_n = cat(1,m(modn,:),m_all_rip_n);
            m_all_rip = cat(1,m,m_all_rip);
            
            
            m = mean_assembly{CSPLUS,1};
            if size(m,2)==1
                m=m';
            end
            
            
            m = bsxfun(@rdivide, m, ass_norm');
%             m_all_stim1= cat(1,m(modhi,:),m_all_stim1);
            m_all_stim1= cat(1,m,m_all_stim1);
            
            m = mean_assembly{CSMINUS,1};
            if size(m,2)==1
                m=m';
            end
            
            m = bsxfun(@rdivide, m, ass_norm');
%             m_all_stim2= cat(1,m(modhi,:),m_all_stim2);
            m_all_stim2= cat(1,m,m_all_stim2);

        end
        
        m_all_rip_area=[];
        for ises_ = ises_area'

            mean_assembly = mean_assembly_all{iarea,ises_,5};

            if isempty(mean_assembly)
                continue
            end
            m = mean_assembly{RIP,1};
            if size(m,2)==1
                m=m';
            end
            m_all_rip_area = cat(1,m,m_all_rip_area);
        
        end
        
        %%
        
        timebins = -1:0.001:1;
        
        idx_pre = find(timebins>=-1 & timebins<=0);
        
        idx = find(timebins>=-0.1 & timebins<=0.1);
        
        % Normalizing stim (zscore)
        norm = cat(2,m_all_stim1,m_all_stim2);
        m = mean(norm(:,[idx_pre idx_pre+size(m_all_stim1,2)]),2);
        s = std(norm(:,[idx_pre idx_pre+size(m_all_stim1,2)]),0,2);

        m_all_stim1 = bsxfun(@rdivide,bsxfun(@minus,m_all_stim1,m),s);
        m_all_stim2 = bsxfun(@rdivide,bsxfun(@minus,m_all_stim2,m),s);
        
        % Normalizing ripples (zscore)
        m = mean(m_all_rip(:,abs(timebins)>0.05),2);
        s = std(m_all_rip(:,abs(timebins)>0.05),0,2);
        m_all_rip = bsxfun(@rdivide,bsxfun(@minus,m_all_rip,m),s);
        
        m = mean(m_all_rip_hi(:,abs(timebins)>0.05),2);
        s = std(m_all_rip_hi(:,abs(timebins)>0.05),0,2);
        m_all_rip_hi = bsxfun(@rdivide,bsxfun(@minus,m_all_rip_hi,m),s);
        
        m = mean(m_all_rip_lo(:,abs(timebins)>0.05),2);
        s = std(m_all_rip_lo(:,abs(timebins)>0.05),0,2);
        m_all_rip_lo = bsxfun(@rdivide,bsxfun(@minus,m_all_rip_lo,m),s);
        
        m = mean(m_all_rip_n(:,abs(timebins)>0.05),2);
        s = std(m_all_rip_n(:,abs(timebins)>0.05),0,2);
        m_all_rip_n = bsxfun(@rdivide,bsxfun(@minus,m_all_rip_n,m),s);
        
        
        %%
        timebins = -1:0.001:4;
        idx_stim = find(timebins>=0 & timebins<=2);
        idx_trace = find(timebins>=2 & timebins<=3);
        if iperiod==1
            idx_period = idx_stim;
        else
            idx_period = idx_trace;
        end
        
        aux = m_all_stim1;

        peak = mean(aux(:,idx_period),2);
        
        idx_rip1 = find(peak>=quantile(peak,0));
        peak1=peak;
        
        stim1= mean(aux(:,idx_stim),2);
        trace1= mean(aux(:,idx_trace),2);
        
        aux = m_all_stim2;
        peak = mean(aux(:,idx_period),2);
        
        idx_rip2 = find(peak>=quantile(peak,0));
        peak2=peak;
        stim2= mean(aux(:,idx_stim),2);
        trace2= mean(aux(:,idx_trace),2);
        
        
        %%
        
        figure(1)
        
        timebins = -1:0.001:1;
        
        aux = m_all_rip_hi(:,:);
%         aux = m_all_rip(idx_rip1,:);
        
        sub(4,3,iperiod+2*(iarea-1),3)
        maux = mean(aux(:,abs(timebins)<0.025),2);
        m=nanmean(maux);
        s=nanstd(maux)./sqrt(sum(~isnan(maux)));
        
        hold on
        bar((ilearn-1)*3+3,m,'facecolor',col(1+3*(ilearn-1),:))
        errorbar((ilearn-1)*3+3,m,s,'.k')
        maux_all{ilearn,1}=maux;
        
        sub(4,3,iperiod+2*(iarea-1),1:2)
        
        hold on


        m = nanmean(aux);
        s = nanstd(aux)./sqrt(sum(~isnan(aux)));
        
        plot(timebins,m,'color',col(1+3*(ilearn-1),:),'linewidth',2)
        hold on
        plot(timebins,m+s,'color',col(1+3*(ilearn-1),:),'linewidth',1)
        plot(timebins,m-s,'color',col(1+3*(ilearn-1),:),'linewidth',1)
        
        
        
        aux = m_all_rip_lo(:,:);
        
        sub(4,3,iperiod+2*(iarea-1),3)
        maux = mean(aux(:,abs(timebins)<0.025),2);
        m=nanmean(maux);
        s=nanstd(maux)./sqrt(sum(~isnan(maux)));
        
        hold on
        bar((ilearn-1)*3+1,m,'facecolor',col(2+3*(ilearn-1),:))
        errorbar((ilearn-1)*3+1,m,s,'.k')
        maux_all{ilearn,3}=maux;
        ylabel('Assembly reactivation')
        
        
        sub(4,3,iperiod+2*(iarea-1),1:2)
        
        m = nanmean(aux);
        s = nanstd(aux)./sqrt(sum(~isnan(aux)));
        
        plot(timebins,m,'color',col(2+3*(ilearn-1),:),'linewidth',2)
        hold on
        plot(timebins,m+s,'color',col(2+3*(ilearn-1),:),'linewidth',1)
        plot(timebins,m-s,'color',col(2+3*(ilearn-1),:),'linewidth',1)
        
        
        aux = m_all_rip_n(:,:);
        
        sub(4,3,iperiod+2*(iarea-1),3)
        maux = mean(aux(:,abs(timebins)<0.025),2);
        m=nanmean(maux);
        s=nanstd(maux)./sqrt(sum(~isnan(maux)));
        
        hold on
        bar((ilearn-1)*3+2,m,'facecolor',col(3+3*(ilearn-1),:))
        errorbar((ilearn-1)*3+2,m,s,'.k')
        maux_all{ilearn,2}=maux;
        
        sub(4,3,iperiod+2*(iarea-1),1:2)
        
        hold on

        m = nanmean(aux);
        s = nanstd(aux)./sqrt(sum(~isnan(aux)));
        
        plot(timebins,m,'color',col(3+3*(ilearn-1),:),'linewidth',2)
        hold on
        plot(timebins,m+s,'color',col(3+3*(ilearn-1),:),'linewidth',1)
        plot(timebins,m-s,'color',col(3+3*(ilearn-1),:),'linewidth',1)
        
        xlim([-0.1 0.1])
        xlabel('Time from RIP (s)')
        ylabel('Assembly Reactivation (+SEM)')
        title([area_labels{iarea} '-' learn_labels{ilearn} ': CS' ])

        if iperiod==1
            title([area_labels{iarea} ' - Modulated Assemblies (STIM)'])
        else
            title([area_labels{iarea} ' - Modulated Assemblies (TRACE)'])
        end
            
        if iarea==1
            ylim([-5 90])
        elseif iarea==2
            ylim([-5 20])
        end
        
        if iarea==1
            xl =[-4 4];
            yl = [-20 150];
            xl2 =[-2.5 2.5];
            yl2 = xl2;
        else
            xl =[-8 8];
            yl = [-5 30];
            xl2 =[-2 14];
            yl2 = [-2 2];
        end
        
        peak_label = {'Stim','Trace'};
        
        %%
        figure(2)
        
        sub(2,4,1,iperiod+(iarea-1)*2)
        hold on 
                
        aux = m_all_rip;
        mrip = mean(aux(:,abs(timebins)<0.025),2);

        invalid = isnan(mrip)|isnan(peak1)|isnan(peak2);
        mrip(invalid) = [];
        peak1(invalid) = [];
        peak2(invalid) = [];
        stim1(invalid) = [];
        stim2(invalid) = [];
        trace1(invalid) = [];
        trace2(invalid) = [];
        
        
        activation1= (trace1-stim1);
        activation2= (trace2-stim2);
        modulation = (peak1-peak2);
        if iperiod==1
            activation=activation1;
        elseif iperiod==2
            activation =activation2;
        end
        
        activation_all = [activation_all; activation];
        
        modulation_all = [modulation_all; modulation];
        mrip_all = [mrip_all; mrip];
        
        [sr1(ilearn) sp1(ilearn)] = corr(modulation,mrip,'type','Spearman');

        plot(modulation,mrip,'.','color',col2(ilearn,:))
        xlim(xl)
        ylim(yl)
        xlabel(['Trial-tipe modulation ' peak_label{iperiod}])
        ylabel('SWR Reactivation')
        
        cs_label = {'CS+', 'CS-'};
        sub(2,4,2,iperiod+(iarea-1)*2)
        plot(activation,mrip,'.','color',col2(ilearn,:))
        hold on 
        xlim(xl)
        ylim(yl)
        xlabel(['Trial-period modulation (Stim-Trace) - ' cs_label{iperiod}])
        ylabel('SWR Reactivation')
        
    end
    
    %%
    figure(2)
        
    sub(2,4,1,iperiod+(iarea-1)*2)
    plot([0 0], ylim,'k--')
    [sr1 sp1] = corr(modulation_all,mrip_all,'type','Spearman');
    title({area_labels{iarea} ['p=' num2str(sp1)], ['R:' num2str(sr1)]})
    
    sub(2,4,2,iperiod+(iarea-1)*2)
    plot([0 0], ylim,'k--')
    [sr1 sp1] = corr(activation_all,mrip_all,'type','Spearman');
    title({area_labels{iarea} ['p=' num2str(sp1)], ['R:' num2str(sr1)]})
   
    
    
    figure(1)
end
end

sub(4,3,4,1:2)
legend({'CS- suppressed pre','CS+ suppressed learning','Non-modulated pre'...
    'CS- suppressed post','CS+ suppressed post','Non-modulated post'})

figure(2)
legend({'Pre' 'Post'})

%%



%% Figure 6D

clear
load('preprocessed_data', 'data', 'hasPF', 'hasCA1', 'learned', 'ses', 'exp', 'trainingday', 'animal',...
    'varnames', 'srate', 'datafiles')

%%

day=[];
animal=[];
for i=1:length(data)
    day(i)= data(i).TrainingDay{1}(1);
    animal(i)= data(i).Animal{1}; 
end

nlick_all={};
for ises = 1:length(data)
    disp(['Performing session: ' num2str(ises)])
    animal(ises) = data(ises).Animal{1};
    spktimes = double(data(ises).Spikes{1})/srate;
    spkid = data(ises).Spikes{2};
    area = data(ises).Area{1};
    etime = data(ises).Events{1}/srate;
    etype = data(ises).Events{2};
    spkid = spkid-min(spkid)+1;
    classes = unique(spkid);

    %%
    idx1 = find(etype == 1);
    trialtype = 1*ones(length(idx1),1);
    idx2 = find(etype == 2);
    trialtype = [trialtype; 2*ones(length(idx2),1)];
    idx = [idx1; idx2];
    trialtime = etime(idx);
    trialtime(:,2) = etime(idx)+3;
    
    valid = 0;

    step=200/1000;
    timebins = -1:step:4;
    binwidth=step;
    preidx = find(timebins>-1 & timebins+binwidth<=0);
    traceidx = find(timebins>2 & timebins+binwidth<=3);
    
    nlick=[];
    nlick_isi=[];
    nlick_first=[];
    for itrial = 1:size(trialtime,1)
        valid = etime>=trialtime(itrial,1)-1 & etime<trialtime(itrial,2)+1;
        valid = valid & etype==3;

        etime_ = etime(valid);
        etime_ = etime_-trialtime(itrial,1);
        etype_ = etype(valid);
        
        for ibin=1:length(timebins)-1
            idx = find(etime_>timebins(ibin) & etime_<=timebins(ibin+1));
            nlick(itrial,ibin) = length(idx);
            
        end
        idx = find(etime_>0 & etime_<=2);
        nlick_isi(itrial,1) = mean(diff(etime_(idx)));
        
        nlick_first(itrial,1)=NaN;
        if ~isempty(idx)
            nlick_first(itrial,1) = min(etime_(idx));
        end
        
        idx = find(etime_>2 & etime_<=3);
        nlick_isi(itrial,2) = mean(diff(etime_(idx)));
        
        nlick_first(itrial,2)=NaN;
        if ~isempty(idx)
            nlick_first(itrial,2) = min(etime_(idx));
        end
        
        idx = find(etime_>3 & etime_<=4);
        nlick_isi(itrial,3) = mean(diff(etime_(idx)));
        nlick_first(itrial,3)=NaN;
        if ~isempty(idx)
            nlick_first(itrial,3) = min(etime_(idx));
        end
        
        
        %%
    end
    trialtype_all{ises} = trialtype;
    trialtime_all{ises} = trialtime;
    nlick_all{ises} = nlick;
    nlick_isi_all{ises} =nlick_isi;
    nlick_first_all{ises} =nlick_first;
end


idx_pre = find(timebins<-0.1);
idx_stim = find(timebins>=-0.1 & timebins<2);
idx_trace = find(timebins>=2 & timebins<3);
idx_rwd = find(timebins>=3);
%%

CA1=1;
PFC=2;
            
PRE = 1;
STIM = 2;
TRACE = 3;
RWD = 4;
ALL = 5;

timebins = -1:0.001:4;
periods={};
periods{PRE} = find(timebins<0);
periods{STIM} = find(timebins>=0 & timebins<1);
periods{TRACE} = find(timebins>=2 & timebins<3);
periods{RWD} = find(timebins>=3 & timebins<4);
periods{ALL} = find(timebins>=-Inf);


doplot = false;
valid_ses=[hasCA1', hasPF', hasCA1' & hasPF'];
valid_neurons={CA1, PFC, []};
nsub = 20;
binwitdh = 20;
nneurons_all = [];
ntrials_all = [];

area_labels={'CA1', 'PFC', 'CA1 & PFC'};

%%
clear mean_assembly_all patterns_all n_rip_all quart_all
%%
mean_assembly_all={};


try
    load(['Assembly_patterns'])
    computeassemblies=0;
catch
    computeassemblies=1;
end

%%
try
    load(['Figure6_data'],...
        'n_rip_all','quart_all')

catch
    clear quart_all mod
    %%
    for ises = find(valid_ses(:,3))'

        %seed = 654299090;
        %rng(seed)

        clear mean_assembly

        animal(ises) = data(ises).Animal{1};
        spktimes = double(data(ises).Spikes{1})/srate;
        spkid = data(ises).Spikes{2};
        area = data(ises).Area{1};
        etime = data(ises).Events{1}/srate;
        etype = data(ises).Events{2};
        spkid = spkid-min(spkid)+1;
        classes = unique(spkid);
        %%
        if isempty(data(ises).RipTS{1})
            continue
        end
        %%
        riptime = double(data(ises).RipTS{1})/srate;

        area_aux =[];
        area_aux(classes) =area;


        lessneuron=[];


        idx1 = find(etype == 1);
        trialtype = 1*ones(length(idx1),1);
        idx2 = find(etype == 2);
        trialtype = [trialtype; 2*ones(length(idx2),1)];
        idx = [idx1; idx2];
        trialtime = etime(idx);
        trialtime(:,2) = etime(idx)+3;

        timebins = -1:0.001:4;    
        %%
        assemblies_CA1 =[];
        assemblies_PFC =[];
        assemblies_rip_CA1=[];

        for iarea = [CA1 PFC]

            disp(area_labels{iarea})

            nsubsample=1;

            for isample = 1:nsubsample
                area_=area;
                area_both=area;

                valid = 0;

                timebins = -1:0.001:4;

                for isource=ALL

                P = patterns_all{iarea,ises,isource};
                %%
                timebins = -1:0.001:4;
                psth=zeros(length(classes),size(trialtime,1),length(timebins));
                psth_smooth=psth;
                for itrial = 1:size(trialtime,1)
                    valid = spktimes>=trialtime(itrial,1)-1 & spktimes<trialtime(itrial,2)+1;

                    if ~isempty(valid_neurons{iarea})
                        area_aux =[];
                        area_aux(classes) =area_;
                        valid = valid & (area_aux(spkid)==valid_neurons{iarea})';
                    else
                        area_aux =[];
                        area_aux(classes) =area_both;
                        valid = valid & (area_aux(spkid)>0)';
                    end

                    spktimes_ = spktimes(valid);
                    spktimes_ = spktimes_-trialtime(itrial,1);

                    spkid_ = spkid(valid);

                    for ineuron=1:length(classes)
                        idx = find(spkid_==classes(ineuron));
                        if isempty(idx)
                            continue
                        end
                        count = histc(spktimes_(idx),timebins);
                        psth(ineuron,itrial,:)=count;
                        psth_smooth(ineuron,itrial,:)=conv(count,rectwin(binwitdh)/sum(binwitdh),'same');

                    end
                end


                %%

                nlick = nlick_all{ises};
                lick_trace = mean(nlick(:,idx_trace),2);

                ntrial = min(sum(trialtype==1),sum(trialtype==2));

                right1 = find(trialtype==1 & lick_trace > 0);
                wrong1 = find(trialtype==1 & lick_trace <= 0);
                wrong2 = find(trialtype==2 & lick_trace > 0);
                right2 = find(trialtype==2 & lick_trace <= 0);

                %%

                psth_proj=zeros(size(P,2),size(psth_smooth,2),length(timebins));
                parfor itrial=1:size(psth_smooth,2)
                    psth_proj(:,itrial,:) = assembly_activity(P,squeeze(psth_smooth(:,itrial,:)));
                end            

                if isempty(P) | size(P,2)==1
                    continue
                end

                psth_proj_temp = psth_proj;
                mean_proj = squeeze(mean(psth_proj,2));
                mean_proj1 = squeeze(mean(psth_proj(:,trialtype==1,:),2));
                mean_proj2 = squeeze(mean(psth_proj(:,trialtype==2,:),2));

                CSmodulation = mean(mean_proj1(:,periods{TRACE}),2)-mean(mean_proj2(:,periods{TRACE}),2);
                STIMmodulation = mean(mean_proj2(:,periods{TRACE}),2)-mean(mean_proj2(:,periods{STIM}),2);

                if iarea==1
                    [v iCA1] = sort(CSmodulation);
                    iCA1 = iCA1([1 1 end]);
                    [v idx] = min(abs(CSmodulation));
                    iCA1(2) = idx;
                    assemblies_CA1 = psth_proj(iCA1,:,:);
                    if iCA1(1)==iCA1(2)
                        assemblies_CA1(1,:,:)=nan;
                    end
                    if iCA1(3)==iCA1(2)
                        assemblies_CA1(3,:,:)=nan;
                    end

                elseif iarea==2

                    [v iPFC] = sort(STIMmodulation);
                    iPFC = iPFC([1 1 end]);
                    [v idx] = min(abs(STIMmodulation));
                    iPFC(2) = idx;
                    assemblies_PFC = psth_proj(iPFC,:,:);
                    if iPFC(1)==iPFC(2)
                        assemblies_PFC(1,:,:)=nan;
                    end
                    if iPFC(3)==iPFC(2)
                        assemblies_PFC(3,:,:)=nan;
                    end
                end



                %%

                timebins = -1:0.001:1;
                riptime = double(data(ises).RipTS{1})/srate;
                riptime=riptime(:);
                riptime(:,2) = riptime+1;
                rip_prevtrial = zeros(1,size(riptime,1));
                rip_nexttrial = zeros(1,size(riptime,1));
                clear psth psth_smooth
                psth = zeros(length(classes),size(riptime,1),length(timebins));
                psth_smooth = zeros(length(classes),size(riptime,1),length(timebins));

                for irip = 1:size(riptime,1)

                    valid = spktimes>=riptime(irip,1)-1 & spktimes<riptime(irip,2);


                    [v1,idx] = sort(trialtime(:,1));
                    [v2,iidx] = sort(idx);

                    aux = find(riptime(irip,1)<=trialtime(idx(2:end),1)-1 & riptime(irip,1)>trialtime(idx(1:end-1),2)+1);

                    if ~isempty(aux)
                        rip_prevtrial(irip) = find(iidx==aux);
                        rip_nexttrial(irip) = find(iidx==aux+1);
                    end


                    if ~isempty(valid_neurons{iarea})
                        area_aux =[];
                        area_aux(classes) =area_;
                        valid = valid & (area_aux(spkid)==valid_neurons{iarea})';
                    else
                        area_aux =[];
                        area_aux(classes) =area_both;
                        valid = valid & (area_aux(spkid)>0)';
                    end
                    spktimes_ = spktimes(valid);
                    spktimes_ = spktimes_-riptime(irip,1);

                    spkid_ = spkid(valid);

                    for ineuron=1:length(classes)
                        idx = find(spkid_==classes(ineuron));
                        if isempty(idx)
                            continue
                        end
                        count = histc(spktimes_(idx),timebins);
                        psth(ineuron,irip,:)=count;
                        psth_smooth(ineuron,irip,:)=conv(count,rectwin(binwitdh),'same');

                    end
                    %%
                end

                psth_proj_rip=zeros(size(P,2),size(psth_smooth,2),length(timebins));
                parfor irip=1:size(psth_smooth,2)
                    psth_proj_rip(:,irip,:) = assembly_activity(P,squeeze(psth_smooth(:,irip,:)));
                end                        
                window = find(timebins>-0.025 & timebins<0.025);

                m_all_rip = squeeze(nanmean(psth_proj_rip(:,:,:),2));

                norm = m_all_rip(:,timebins<=-0.05 | timebins>=0.05);
                z_m_all_rip  = bsxfun(@rdivide, bsxfun(@minus, m_all_rip, mean(norm,2)), std(norm,0,2));


                if iarea==1
                    ripples_CA1 = squeeze(nanmean(psth_proj_rip(:,:,window),3));
                    mwin = mean(z_m_all_rip(:,window),2);
                    [v iCA1_rip] = sort(mwin);

                    iCA1_rip = iCA1_rip([1 1 end]);
                    [v idx] = min(abs(mwin));
                    iCA1_rip(2) = idx;
                    assemblies_rip_CA1 = psth_proj_temp(iCA1_rip,:,:);
                elseif iarea==2
                    ripples_PFC = squeeze(nanmean(psth_proj_rip(:,:,window),3));
                end

                end
            end
        end

        %%
        n_rip_all(ises) = length(ripples_CA1);

        if length(assemblies_CA1)==0 | length(assemblies_PFC)==0
            continue
        end
        for i=1:size(ripples_CA1,1)
            for j=1:size(ripples_PFC,1)

                [R, p] = corr(ripples_CA1(i,:)',ripples_PFC(j,:)');

                R_rip_all(ises,i,j) = R;
                p_rip_all(ises,i,j) = p;

                t1=quantile(ripples_CA1(i,:),0.5);
                t2=quantile(ripples_PFC(j,:),0.5);

                quart = sum(ripples_CA1(i,:)>t1 & ripples_PFC(j,:)>t2)/length(ripples_CA1);
                quart_all(ises,i,j,1) = quart;

                quart = sum(ripples_CA1(i,:)<=t1 & ripples_PFC(j,:)<=t2)/length(ripples_CA1);
                quart_all(ises,i,j,2) = quart;

                quart = sum(ripples_CA1(i,:)>t1 & ripples_PFC(j,:)<=t2)/length(ripples_CA1);
                quart_all(ises,i,j,3) = quart;

                quart = sum(ripples_CA1(i,:)<=t1 & ripples_PFC(j,:)>t2)/length(ripples_CA1);
                quart_all(ises,i,j,4) = quart;

            end
        end

    end
    save(['Figure6_data'],...
        'n_rip_all','quart_all','-append')

end

%%


figure(3);clf

xl=[0 5];
yl=[0 .5];

col=[.2 .2 .2;
    .2 .7 .2];

for ilearn=1:2
    ises = find(valid_ses(:,3) & learned(:,1)==ilearn-1);
    
    quart = quart_all(ises,:,:,:);
    quart(quart==0)= nan;
    quart = permute(quart,[4 1 2 3]);
    
    quart = quart(:,:);
    hold on
    plot([0 0],[0 0],'.','color',col(1,:))
    plot([0 0],[0 0],'.','color',col(2,:))

    for i=[1 3]
        plot((ilearn-1)*3+i/2-0.1+0.1*randn(1,length(quart(i,:))),quart(i,:),'.','color',col(ilearn,:))
    end
    
    ylabel('% of ripples')
    xlabel('Quadrant')
    set(gca,'xtick',[0.5 1.5 3.5 4.5],'xticklabel',{'Q1','Q2'})
    ylim(yl)
    xlim(xl)
    
end

plot(xlim,0.25*[1 1],'k--')
legend({'Pre' 'Post'})


