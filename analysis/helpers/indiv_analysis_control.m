function indiv_analysis_control(subj,refit)

%% define folders
datafiles = dir(fullfile(datafolder,sprintf('%s*_data.mat',subj)));

%% load data
datas = {};
for fi = 1:length(datafiles)
    load(fullfile(datafolder,datafiles(fi).name));
    % we only want to keep fixation task data, so check that task == 0
    if all(data.design(:,9)==0)
        datas{end+1}=data;
    end
end
    
%% Fit model to all available data

if refit 
    %% fit full model HRF
    fithrf = fitCCTimecourseROIModel(datas,'fithrf');
%     save(fullfile(datafolder,sprintf('%s_fithrf.mat',subj)),'fithrf'); close all

    %%
    deconCohCon(datas,fithrf,'lV1',subj);
    %% fit full model ROIs
    load(fullfile(datafolder,sprintf('%s_fithrf.mat',subj)));
    %%
    fitroi = fitCCTimecourseROIModel(datas,'fitroi',fithrf); close all
    
    %%
    deconCohCon(datas,fitroi,'V1',subj);

    %%
    save(fullfile(datafolder,sprintf('%s_fitroi.mat',subj)),'fitroi');
    close all
    

    %% complete refit (only HRF, but for all ROIs)
%     disp('Running refit'); % this has VERY LITTLE effect... probably just over-fitting
    fitroi = fitCCTimecourseROIModel(datas,'refit',fitroi); close all
%     save(fullfile(datafolder,sprintf('%s_fitroi.mat',subj)),'fitroi');
end

%% fitCCBehavModelROI
% this code estimates a sigma for each ROI combination when fit to the
% behavioral data
% % % % % % % % subjb = subj([1 3 4 5]);
% % % % % % % % adata = loadadata(subjb);
rois = {'V1','V2','V3','V4','V3a','V3b','V7','MT'};
% % % % % % % % asigma = zeros(8,8);
% % % % % % % % likes = zeros(8,8);
% % % % % % % % parfor ri = 1:length(rois)
% % % % % % % %     values = asigma(ri,:);
% % % % % % % %     lvals = likes(ri,:);
% % % % % % % %     for ji = 1:length(rois)
% % % % % % % %         fit_s = fitCCBehavModelROI_control(adata,0,subj,subjb,rois{ri},rois{ji});
% % % % % % % %         values(ji) = fit_s.params.sigma;
% % % % % % % %         lvals(ji) = fit_s.likelihood;
% % % % % % % %     end
% % % % % % % %     asigma(ri,:) = values;
% % % % % % % %     likes(ri,:) = lvals;
% % % % % % % % end
% % % % % % % % fname = fullfile(datafolder,'sigmas.mat');
% % % % % % % % if exist(fname,'file')==2, load(fname); end
% % % % % % % % sigmas.(subj) = asigma;
% % % % % % % % save(fname,'sigmas');
% 
%% deconCohCon plot
load(fullfile(datafolder,sprintf('%s_fitroi.mat',subj)));
for ri = 1:length(rois)
    deconCohCon(datas,fitroi,rois{ri},subj);
    close all
end
% 

%% Output response and model
% f = figure;
load(fullfile(datafolder,sprintf('%s_fitroi.mat',subj)));
load(fullfile(datafolder,sprintf('%s_decon.mat',subj)));
% we're going to plot the model fits onto the actual data for each of these
% rois (we'll do it one plot, because... I don't know whatever)

% cmap = brewermap(7,'PuOr');
for ri = 1:length(rois)
%     subplot(length(rois),1,ri), hold on
    % we'll plot the mean, technically this is incorrect (since it averages
    % across the other feature), but fuck it, we only really care about the
    % shape matching
    ucon = unique(decondata.(rois{ri}).cc.conidxs);
    ucoh = unique(decondata.(rois{ri}).cc.cohidxs);
    for ui = 1:length(ucon)
        resp = decondata.(rois{ri}).cc.resp(logical((decondata.(rois{ri}).cc.cohidxs==0).*(decondata.(rois{ri}).cc.conidxs==ucon(ui))));
        mresp = decondata.(rois{ri}).cc.mresp(logical((decondata.(rois{ri}).cc.cohidxs==0).*(decondata.(rois{ri}).cc.conidxs==ucon(ui))));
%         plot(ucon(ui)*100,mean(resp),'o','MarkerFaceColor',cmap(2,:),'MarkerEdgeColor',[1 1 1],'MarkerSize',10);
        decondata.(rois{ri}).rcon(ui) = resp;
        decondata.(rois{ri}).rconm(ui) = mresp;
    end
    % we also save this subject's responses for this ROI, we'll use these
    % to make averages later
    decondata.(rois{ri}).ucon = ucon;
    for ui = 1:length(ucoh)
        resp = decondata.(rois{ri}).cc.resp(logical((decondata.(rois{ri}).cc.conidxs==0.25).*(decondata.(rois{ri}).cc.cohidxs==ucoh(ui))));
        mresp = decondata.(rois{ri}).cc.mresp(logical((decondata.(rois{ri}).cc.conidxs==0.25).*(decondata.(rois{ri}).cc.cohidxs==ucoh(ui))));
%         plot(ucoh(ui)*100,mean(resp),'o','MarkerFaceColor',cmap(6,:),'MarkerEdgeColor',[1 1 1],'MarkerSize',10);
        decondata.(rois{ri}).rcoh(ui) = resp;
        decondata.(rois{ri}).rcohm(ui) = mresp;
    end
    decondata.(rois{ri}).ucoh = ucoh;
%     x = 0:.01:1;
%     roinums = cellfun(@(x) ~isempty(strfind(x,rois{ri})),fitroi.ROIs,'UniformOutput',false);
%     roinums = find([roinums{:}]);
%     con = mean([conModel(x,fitroi.roiparams{roinums(1)}); conModel(x,fitroi.roiparams{roinums(2)})]);
%     con = con-con(find(x==0.25,1));
%     coh = mean([cohModel(x,fitroi.roiparams{roinums(1)}); cohModel(x,fitroi.roiparams{roinums(2)})]);
%     coh = coh-coh(1);
%     plot(x(x>=0.25)*100,con(x>=0.25),'-','Color',cmap(2,:));
%     plot(x*100,coh,'-','Color',cmap(6,:));
% %     axis([0 1 -0.05 0.2])
%     xlabel('Stimulus strength (%)');
%     ylabel('Signal change (%)');
%     title(rois{ri})
%     drawPublishAxis
end
save(fullfile(datafolder,sprintf('%s_decon.mat',subj)),'decondata');
% 
% savepdf(f,fullfile(datafolder,sprintf('%s_roifmri.pdf',subj)));