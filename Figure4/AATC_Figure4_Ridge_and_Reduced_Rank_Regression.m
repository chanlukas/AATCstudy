% % Code for reproducing Figure 4 and associated Supp. Figures
% % More information on the method used here can be found on Semedo et al., 2019.
% % Uses sub.m auxiliary function to plot.

clear
load('preprocessed_data', 'data', 'hasPF', 'hasCA1', 'learned', 'ses', 'exp', 'trainingday', 'animal',...
    'varnames', 'srate', 'datafiles')

%%

CA1=1;
PFC=2;

PRE=1;
STIM=2;
TRACE=3;
RWD=4;

valid_ses=[hasCA1', hasPF', hasCA1' & hasPF'];
valid_neurons={CA1, PFC, []};

nsub = 20;
area_labels={'CA1', 'PFC', 'CA1 & PFC'};

n_predictive_dim=[];
max_performance=[];
max_performance_ridge=[];

%%
try
    load('Figure4_data','n_predictive_dim','max_performance','nfold',...
        'nsub','lambvec','max_lamb_ridge','max_performance_ridge')
catch
    for ises = find(valid_ses(:,3)==1)'
        %%% Analyses on the paper were done using this random seed:
        % seed = 654299090;
        % rng(seed)

        disp(['Performing session: ' num2str(ises)])

        animal(ises) = data(ises).Animal{1};
        spktimes = double(data(ises).Spikes{1})/srate;
        spkid = data(ises).Spikes{2};
        area = data(ises).Area{1};
        etime = data(ises).Events{1}/srate;
        etype = data(ises).Events{2};
        spkid = spkid-min(spkid)+1;
        classes = unique(spkid);

        area_aux =[];
        area_aux(classes) =area;


        nneurons=[];
        for iarea=1:length(area_labels)-1
            nneurons(iarea) = sum(area_aux(unique(spkid))==valid_neurons{iarea});
        end

        if sum(nneurons==0)
            nneurons(iarea+1) = 0;
            lessneuron = [];
        else
            nneurons(iarea+1)= sum(nneurons(1:2));
            lessneuron = min(nneurons);
            if mod(lessneuron,2)==1
                lessneuron=lessneuron-1;
            end
        end

        %% Finding CS+ and CS- times

        idx1 = find(etype == 1);
        trialtype = 1*ones(length(idx1),1);
        idx2 = find(etype == 2);
        trialtype = [trialtype; 2*ones(length(idx2),1)];
        idx = [idx1; idx2];
        trialtime = etime(idx);
        trialtime(:,2) = etime(idx)+3;

        timebins = -1:0.001:4;


        %%
        psth_rrr = {};

        for iarea = 1:2 % 1 = CA1; 2 = PFC

            disp(area_labels{iarea})

            %%

            psth=zeros(length(classes),size(trialtime,1),length(timebins));
            valid = 0;


            %%

            for itrial = 1:size(trialtime,1)
                valid = spktimes>=trialtime(itrial,1)-1 & spktimes<trialtime(itrial,2)+1;
                valid = valid & (area_aux(spkid)==valid_neurons{iarea})';

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

                end
                %%
            end

            %%
            binwidth=200;
            step = 200;
            newbins = 1:step:length(timebins)-binwidth;
            psth_aux = [];
            for ibin = 1:length(newbins)
                psth_aux(:,:,ibin) = sum(psth(:,:,newbins(ibin):newbins(ibin)+binwidth),3);
            end
            psth_rrr{iarea} = psth_aux(find(sum(sum(psth,3),2)~=0),:,:);
            %%
        end

        %%    


        % Computing Residuals
        X_all = bsxfun(@minus,psth_rrr{1}, mean(psth_rrr{1},2));
        Y_all = bsxfun(@minus,psth_rrr{2}, mean(psth_rrr{2},2));

        % Computing training and test set sizes. Equal sizes for both areas 
        maxn = min(size(X_all,1),size(Y_all,1));
        ntrain = round(0.5*maxn);
        ntest = maxn-ntrain;


        bincol=1:4;
        col=[.7 .7 .7;
             .7 .7 .2;
             .2 .2 .7;
             .2 .7 .2];
        if size(X_all,1)<=ntrain | size(Y_all,1)<=ntrain 
            continue
        end

        nsub = 20; % number of subsamples (without replacement)
        for isub=1:nsub

            % Defining Target and Source population from each area
            [v,idx] = sort(rand(1,size(X_all,1)));
            X_source = X_all(idx(1:ntrain),:,:);
            X_target = X_all(idx((1:ntest)+ntrain),:,:);

            [v,idx] = sort(rand(1,size(Y_all,1)));
            Y_source= Y_all(idx(1:ntrain),:,:);
            Y_target= Y_all(idx((1:ntest)+ntrain),:,:);

            sources = {X_source, Y_source}; % X = CA1; Y = PFC
            targets = {X_target, Y_target};

            % Computing Reduced Rank Regression and Ridge-Regression
            % with the different source-target pairs

            bingroup = {1:5, 6:15, 16:20,21:25}; % psth indices for each trial period (PRE, STIM, TRACE and RWD)
            for isource = 1:length(sources)
                for itarget= 1:length(targets)

                    clear m_r s_r m_rrr s_rrr 
                    for ibin=[PRE, STIM, TRACE, RWD]
                        X=squeeze(sources{isource}(:,:,bingroup{ibin}));
                        Y=squeeze(targets{itarget}(:,:,bingroup{ibin}));

                        id = repmat(trialtype,1,size(X,3));
                        id = reshape(id,1,[]);

                        clear performance

                        % Computing residuals
                        X = reshape(X,size(X,1),[])';
                        Y = reshape(Y,size(Y,1),[])';

                        X = bsxfun(@minus,X,mean(X,1));
                        Y = bsxfun(@minus,Y,mean(Y,1));


                        n=size(X,1);
                        p = size(X,2);
                        q = size(Y,2);
                        nfold = 10;

                        idx = crossvalind('Kfold',n,nfold);
                        ridge_performance=[];

                        % 10-fold crossvalidation
                        for ifold = 1:nfold
                            X_train = X(idx~=ifold,:);
                            X_test = X(idx==ifold,:);
                            Y_train = Y(idx~=ifold,:);
                            Y_test = Y(idx==ifold,:);

                            id_train = id(idx~=ifold);
                            id_test = id(idx==ifold);


                            % In the model defined by Y=XB, let's find B using the ordinary least 
                            % square solution (Bols)

                            % Ridge regression
                            lambvec = 0:0.01:5; 
                            for ilamb = 1:length(lambvec)
                                lambda = lambvec(ilamb);
                                Bols = pinv((X_train'*X_train)+lambda*n*eye(p))*X_train'*Y_train;    
                                aux = X_train*Bols;

                                [coeff score] = pca(aux);

                                V = coeff(:,:); % use all the ranks

                                B_ = Bols*V;
                                Y_ = X_test*B_*V'; %Yrrr

                                SSres = sum(reshape((Y_test-Y_).^2,1,[]));
                                SStot = sum(reshape(bsxfun(@minus,Y_test,mean(Y_test)).^2,1,[]));
                                R2 = 1-((SSres)/(SStot));

                                ridge_performance(ifold,ilamb) = R2;
                            end

                            % Reduced rank regression
                            rank = [1:10 size(coeff,2)];
                            for irank = 1:length(rank)

                                if rank(irank)>size(coeff,2)
                                    RRRperformance(ifold,irank) = nan;
                                    continue                                
                                end

                                lambda=0; % no regularization

                                Bols = pinv((X_train'*X_train)+lambda*n*eye(p))*X_train'*Y_train;

                                V = coeff(:,1:rank(irank)); % reducing the ranks

                                B_ = Bols*V;
                                Y_ = X_test*B_*V'; %Yrrr


                                SSres = sum(reshape((Y_test-Y_).^2,1,[]));
                                SStot = sum(reshape(bsxfun(@minus,Y_test,mean(Y_test)).^2,1,[]));
                                R2 = 1- ((SSres)/(SStot));

                                RRRperformance(ifold,irank) = R2;
                            end

                        end

                        % Finding the optimal lambda for the Ridge regression
                        mrid = mean(ridge_performance);
                        srid = std(ridge_performance)./sqrt(sum(~isnan(ridge_performance)));

                        [v, imax] = max(mrid);
                        lamb = find(mrid>v-srid(imax),1);

                        max_lamb_ridge(ises,isub,isource,itarget,ibin) = lamb;
                        max_performance_ridge(ises,isub,isource,itarget,ibin) = mrid(lamb);

                        %%

                        m_r(ibin,:) = mrid;
                        s_r(ibin,:) = srid;

                        m_rrr(ibin,:) = mean(RRRperformance);
                        s_rrr(ibin,:) = std(RRRperformance)./sqrt(sum(~isnan(RRRperformance)));
    %                   

                    end

                    % Extracting the predictive dim. and the maximal
                    % performance for RRR.
                    for ibin=1:size(m_r,1)
                        [v, imax] = max(m_rrr(ibin,:));
                        n_predictive_dim(ises,isub,isource,itarget,ibin) = find(m_rrr(ibin,:)>v-s_rrr(ibin,imax),1);
                        max_performance(ises,isub,isource,itarget,ibin) = m_rrr(ibin,n_predictive_dim(ises,isub,isource,itarget,ibin));
                    end
                end
            end

        end     
    end

    %%
    save('Figure4_data','n_predictive_dim','max_performance','nfold',...
        'nsub','lambvec','max_lamb_ridge','max_performance_ridge')
end


%% Ploting performance of Ridge Regression

figure(1); clf
col = [.7 .7 .7;
    .2 .7 .2];

nbin=4;
for isource=1:2
    for itarget=1:2
        for ilearn=1:2
            ises=find(valid_ses(:,3) & learned(:,1)==ilearn-1);    
            clear m s 
            aux = squeeze(mean(max_performance_ridge(ises,:,isource,itarget,:),2));

            m = mean(aux,1);
            s = std(aux,0,1)./sqrt(sum(~isnan(aux)));
            sub(2,2,1,isource)
            hold on
            
            if ilearn==1
                bar((1:nbin)*2-2+ilearn+(itarget-1)*10,m,0.4,'facecolor', col(1,:))
                bar((1:nbin)*2-2+ilearn+(itarget-1)*10,m,0.4,'facecolor', col(2,:))
            end
            bar((1:nbin)*2-2+ilearn+(itarget-1)*10,m,0.4,'facecolor', col(ilearn,:))
            errorbar((1:nbin)*2-2+ilearn+(itarget-1)*10,m,s,'.k')
            set(gca,'xtick', [1.5:2:9 11.5:2:19], 'xticklabel',{'pre', 'stim', 'trace','rwd'},'XTickLabelRotation',45)
            xlabel({'Target:'; 'CA1        /        PFC'})
            ylabel('R2')
        end
        
    end
    
    
    sub(2,2,1,isource)
    legend({'pre' 'post'})
    ylim([0 5])
    ylim([0 0.4])
    title({['Source activity: ' area_labels{isource}]})
    
end

%% Ploting predictive dim. and performance of RRR

figure(2); clf
col = [.7 .7 .7;
    .2 .7 .2];

nbin=4;
for isource=1:2
    p=[];
    p2=[];
    paux={};
    for itarget=1:2
        for ilearn=1:2
            ises=find(valid_ses(:,3) & learned(:,1)==ilearn-1);    
            
            aux = squeeze(mean(max_performance(ises,:,isource,itarget,:),2));

            m = mean(aux,1);
            s = std(aux,0,1)./sqrt(sum(~isnan(aux)));
            
            sub(2,2,1,isource)
            hold on
            if ilearn==1
                bar((1:nbin)*2-2+ilearn+(itarget-1)*10,m,0.4,'facecolor', col(1,:))
                bar((1:nbin)*2-2+ilearn+(itarget-1)*10,m,0.4,'facecolor', col(2,:))
            end
            bar((1:nbin)*2-2+ilearn+(itarget-1)*10,m,0.4,'facecolor', col(ilearn,:))
            errorbar((1:nbin)*2-2+ilearn+(itarget-1)*10,m,s,'.k')
            set(gca,'xtick', [1.5:2:9 11.5:2:19], 'xticklabel',{'pre', 'stim', 'trace','rwd'},'XTickLabelRotation',45)
            xlabel({'Target:'; 'CA1        /        PFC'})
            ylabel('R2')
            
            %%
            aux = squeeze(mean(n_predictive_dim(ises,:,isource,itarget,:),2));

            m = mean(aux,1);
            s = std(aux,0,1)./sqrt(sum(~isnan(aux)));
            
            sub(2,2,2,isource)
            hold on
            
            if ilearn==1
                bar((1:nbin)*2-2+ilearn+(itarget-1)*10,m,0.4,'facecolor', col(1,:))
                bar((1:nbin)*2-2+ilearn+(itarget-1)*10,m,0.4,'facecolor', col(2,:))
            end
            bar((1:nbin)*2-2+ilearn+(itarget-1)*10,m,0.4,'facecolor', col(ilearn,:))
            errorbar((1:nbin)*2-2+ilearn+(itarget-1)*10,m,s,'.k')
            set(gca,'xtick', [1.5:2:9 11.5:2:19], 'xticklabel',{'pre', 'stim', 'trace','rwd'},'XTickLabelRotation',45)
            xlabel({'Target:'; 'CA1        /        PFC'})
            ylabel('Predictive Dimension')
            
        end
        
    end
    
    sub(2,2,1,isource)
    legend({'pre' 'post'})
    ylim([0 5])
    ylim([0 0.4])
    title({['Source activity: ' area_labels{isource}]})
    
    sub(2,2,2,isource)
    legend({'pre' 'post'})
    ylim([0 5])
    title({['Source activity: ' area_labels{isource}]})    
    
end








