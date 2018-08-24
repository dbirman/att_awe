%%
%sids = {'s043','s300','s305','s329','s025'};
nSIDs = [305 329 43 25 300 346 343 344 338 340 348];
bSIDs = [345 330 337 335 349 354 353 356 334 352]; % behavior only participants
aSIDs = [nSIDs bSIDs];
sids = {};
for bi = 1:length(aSIDs)
    sids{end+1} = sprintf('s%03.0f',aSIDs(bi));
end
rois = {'V1','V2','V3','V4','V3a','V3b','V7','MT'};

%%
%%%%%%%%%%%%%%%%%%%%% NEURAL %%%%%%%%%%%%%%%%%%%%%%%%%
nSIDs = [305 329 43 25 300 346 343 338 344 348];
% NO ATT DATA FROM 340


%    1        2    3     4     5      6   7    8      9
% stimvol basecon lcon rcon basecoh lcoh rcoh timing task


%% Delete old files
% % warning('Are you sure?');
% % keyboard
% % for ni = 1:length(nSIDs)
% %     files = dir(fullfile(datafolder,sprintf('s%04.0f_a*.mat',nSIDs(ni))));
% %     for fi = 1:length(files)
% %         delete(fullfile(datafolder,files(fi).name));
% %     end
% % end
%%

for ni = 1:length(nSIDs)
    indiv_analysis(sprintf('s%04.0f',nSIDs(ni)),'decon');
end

%% Cross analysis
% Cross analysis deconvolves the 16 combinations of contrast and coherence
% rather than separately deconvolving the contrast and coherence response.
% We rely on the rounding of the increments to do this, otherwise there's
% not really enough trials (<~25)
avg_deconResponseAtt_cross;
fit_deconResponseAtt_cross;
indiv_fitbehav_att_cross;

%% Check the gain and increment parameters from fit_deconResponseAtt_cross

load(fullfile(datafolder,'avg_att_cross_fits.mat'));

for ni = 1:length(nSIDs)
    for ri = 1:8
        offset(ni,ri,1) = attfits{ni}{1}.roifit{ri}.params.(sprintf('%soffset_coh',rois{ri}));
        offset(ni,ri,2) = attfits{ni}{1}.roifit{ri}.params.(sprintf('%soffset_con',rois{ri}));
    end
end

%% Bootstrap the difference
diff_offset = offset(:,:,2)-offset(:,:,1);
ci = bootci(1000,@mean,diff_offset);

for ri = 1:8
    disp(sprintf('%s %0.2f\\%%, 95\\%% CI [%0.2f, %0.2f]; ',rois{ri},mean(diff_offset(:,ri)),ci(1,ri),ci(2,ri)));
end
%% Sensitivity space to readout space analysis
sense_readout;


%% don't know what any of this stuff does


% 
% %% TODO: check bheavior equivalent to non-scanner behavior
% cc_rightchoicecontrol;
% % parfor si = 1:length(sids)
% %     subj_behav_analysis_att(sids{si},'refit');
% %     close all
% % end
% % %%
% % for si = 1:length(sids)
% %     subj_behav_analysis_att(sids{si},'disp');
% %     close all
% % end
% % %%
% % for si = 1:length(sids)
% %     subj_behav_analysis_att(sids{si},'right');
% %     close all
% % end
% 
% %% plot individual response plots
% cmap = brewermap(7,'PuOr');
% x = 0:.01:1;
% for si = 1:length(sids)
%     load(fullfile(datafolder,sprintf('%s_data.mat',sids{si})));
%     ycon = conModel(x,data.fits{1}.params);
%     ycoh = cohModel(x,data.fits{1}.params);
%     h = figure; hold on
%     plot(x,ycon,'-','Color',cmap(2,:));
%     plot(x,ycoh,'-','Color',cmap(6,:));
%     drawPublishAxis
%     savepdf(h,fullfile(datafolder,sprintf('%s_response.pdf',sids{si})));
% end
% close all
% 
% %% plot average con/coh functions and mean function
% avg_behavmodel(sids,1);
% avg_behavmodel(sids,2);
% avg_comparebehavmodels(sids);
% %% disp average parameters
% aparams = struct;
% for si = 1:length(sids)
%     subj = sids{si};
%     load(fullfile(datafolder,sprintf('%s_data.mat',subj)));
%     cfit = data.fits{1}; % only use the con-exp coh-exp for now
%     pfields = fields(cfit.params);
%     for pi = 1:length(pfields)
%         if ~isfield(aparams,pfields{pi}),aparams.(pfields{pi}) = nan(size(sids)); end
%         aparams.(pfields{pi})(si) = cfit.params.(pfields{pi});
%     end
% end
% pfields = fields(aparams);
% for pi = 1:length(pfields)
%     ci = bootci(10000,@nanmean,aparams.(pfields{pi}));
%     aparams.(sprintf('%s_ci',pfields{pi})) = ci;
%     aparams.(pfields{pi}) = mean(aparams.(sprintf('%s_ci',pfields{pi})));
% %     aparams.(sprintf('%s',pfields{pi})) = nanmean(aparams.(pfields{pi}));
% end
% afit.params = aparams;
% cc_rightchoicecontrol(zeros(1,12),afit);
% 
% %% BIC Plot
% BICs = zeros(length(sids),3);
% for si = 1:length(sids)
%     subj = sids{si};
%     load(fullfile(datafolder,sprintf('%s_data.mat',subj)));
%     BICs(si,:) = data.BICs;
% end
% BICs = BICs - repmat(BICs(:,1),1,3);
% 
% nakalin_ci = bootci(1000,@mean,BICs(:,2));
% switchstay_ci = bootci(1000,@mean,BICs(:,3));
% 
% disp(sprintf('On average Naka-Rushton/Linear model suffers a BIC difference of %2.3f to %2.3f',nakalin_ci(1),nakalin_ci(2)));
% %% dispInfo (from cohcon;
% 
% % combine data
% control = zeros(length(sids),2,4);
% for si = 1:length(sids)
%     load(fullfile(datafolder,sprintf('%s_data.mat',sids{si})));
%     control(si,:,:) = data.control;
% end
% 
% control(control<=0) = NaN;
% 
% %% Load number of CONTROL trials
% for si = 1:length(sids)
%     [~,n(si)] = loadadata(sids{si});
% end
% 
% %% Fit group model (to estimate what the true underlying functions are
% % so we can compare to neural data and use as a base for indiv models
% % findbestgroupmodel(control,attend); % we won't use the unattended data for this
% 
% %% Plot group performance plot
% avg_performancecontrol;
% 
% %% Threshold Control Plots
% for si = 1:length(sids)
%     cc_thresholdcontrol(sids{si});
%     close all
% end
% 
% %% Plot (threshold across subjects)
% h = figure, hold on
% cmap = brewermap(7,'PuOr');
% hs = zeros(1,4);
% for si = 1:length(sids)
%     subj = sids{si};
%     
%     load(fullfile(datafolder,sprintf('%s_data.mat',subj)));
%     
%     hs(si) = plot(data.contx, data.cont,'-','Color',cmap(2,:));
%     plot(data.cohtx,data.coht,'-','Color',cmap(6,:));
%     plot([0.325, 0.4 0.55, 0.85], data.control(2,:),'s','MarkerFaceColor',cmap(2,:),'MarkerEdgeColor',[1 1 1],'MarkerSize',10);
%     plot([0.15 0.3 0.45 0.6], data.control(1,:),'o','MarkerFaceColor',cmap(6,:),'MarkerEdgeColor',[1 1 1],'MarkerSize',10);
% end
% set(gca,'XTick',[0 0.5 1],'XTickLabels',[0 50 100]);
% set(gca,'YTick',[0 0.2 0.4],'XTickLabels',[0 20 40]);
% xlabel('Stimulus Strength (%)');
% ylabel('Just Noticeable Difference (%)');
% axis([0 1 0 1.1])
% drawPublishAxis
% 
% % savepdf(h,'avg_thresholdmodel.pdf');
% 
% %% Contrast / Coherence parameter averages (overlay functions + average function)
% 
% %% BIC Comparison Plots
% 
% %% Contrast "leak" plots
% 
% % just grab the control condition parameters and compare to attended
% % conditions
% % b_con_... b_att_...
% values = zeros(length(sids),4);
% pfields = {'beta_control_con_cohw','beta_control_coh_conw','bias'};
% differences = zeros(4,2);
% for si = 1:length(sids)
%     load(fullfile(datafolder,sprintf('%s_data.mat',sids{si})));
%     for i = 1:length(pfields)
%         values(si,i) = data.fits{1}.params.(pfields{i});
%     end
% end
% %%u
% coh_leak = values(:,1);
% con_leak = values(:,2);
% bias = values(:,3);
% 
% con_ci = bootci(10000,@mean,con_leak);
% coh_ci = bootci(10000,@mean,coh_leak);
% bias_ci = bootci(10000,@mean,bias);
% 
% %% run threshold controls
% for si = 1:sids
%     cc_thresholdcontrol(sids{si});
%     close all
% end
% 
% %% optimal behavioral functions
% % two subplots: left is 
% x = 0:.01:1;
% 
% load(fullfile(datafolder,'avg_linked.mat'));
% cmap = brewermap(7,'PuOr');
% 
% h = figure;
% group = {'exp','linear'};
% noise = {'additive','poisson'};
% 
% for gi = 1:2
%     for ni = 1:2
%         subplot(2,2,(gi-1)*2+ni); hold on
%         plot(x*100,mean(bootci(1000,@mean,behav.(group{gi}).(noise{ni}).con)),'Color',cmap(2,:));
%         plot(x*100,behav.(group{gi}).(noise{ni}).coh_,'Color',cmap(6,:));
%         xlabel('Stimulus strength (%)');
%         ylabel('Model response (a.u.)');
%         if ni==1
%             axis([0 100 0 50]);
%             set(gca,'XTick',[0 100],'XTickLabel',[0 100]);
%             set(gca,'YTick',[0 50]);
%         else
%             axis([0 100 0 150]);
%             set(gca,'XTick',[0 100],'XTickLabel',[0 100]);
%             set(gca,'YTick',[0 50 100 150]);
%         end
%         title(sprintf('%s: %s',group{gi},noise{ni}));
%         drawPublishAxis
%     end
% end
% 
% savepdf(h,fullfile(datafolder,'avg_optimalbehav.pdf'));
% 
% %% Table of behav info (sigmas and residuals)
% h = figure;
% set(f,'Position',[300 200 300 150]);
% dat = cell(4,3);
% group = {'exp','linear'};
% noise = {'additive','poisson'};
% jg = {'scale','conres','cohres'};
% % start by converting res to BIC
% bic = {'scale','conresBIC','cohresBIC'};
% 
% params = [5 7]; % noise parameters are the same in both (just sigma)
% for gi = 1:2
%     for ni = 1:2
%         for ci = 1:2
%             RSS = behav.(group{gi}).(noise{ni}).(jg{ci+1});
%             n = 21;
%             behav.(group{gi}).(noise{ni}).(bic{ci+1}) = n*log(RSS/n) + params(gi)*log(n);
%         end
%     end
% end
% 
% for gi = 1:2
%     for ni = 1:2
%         for j = 1:3
%             if j==1
%                 dat{(gi-1)*2+ni,j} = behav.(group{gi}).(noise{ni}).(bic{j});
%             else
%                 dat{(gi-1)*2+ni,j} = behav.(group{gi}).(noise{ni}).(bic{j})-behav.exp.additive.(bic{j});
%             end
%         end
%     end
% end
% columnname =   {'Sigma', 'Contrast BIC - Exp/Additive', 'Coherence BIC - Exp/Additive'};
% columnformat = {'numeric', 'numeric', 'numeric'}; 
% t = uitable('Units','normalized','Position',...
%             [0.05 0.05 0.755 0.87], 'Data', dat,... 
%             'ColumnName', columnname,...
%             'ColumnFormat', columnformat,...
%             'RowName',[]); 
% set(gca,'XTick',[]);
% set(gca,'YTick',[]);
% savepdf(h,fullfile(datafolder,'avg_table.pdf'));
