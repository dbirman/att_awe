%% Response correlation analysis
% load subject data
nSIDs = [305 329 43 25 300 346 343 344 338 340 348];
bSIDs = [345 330 337 335 349 354 353 356 334 352]; % behavior only participants
aSIDs = [nSIDs bSIDs];

restructure_afits;

cdir = pwd;

ROIs = {'V1','V2','V3','V4','V3a','V3b','V7','MT'};

%% Compute readout response function -> convert to fMRI
aouts = {};

for ni = 1:length(nSIDs)
    % load subject model fit
    nfit = afits{find(aSIDs==nSIDs(ni))};
    nfit = nfit{1}; % drop poisson model
    % load 
    files = dir(fullfile(datafolder,sprintf('s%04.0f*_data.mat',nSIDs(ni))));
    datas = {};
    for fi = 1:length(files)
        loaded = load(fullfile(datafolder,files(fi).name));
        datas{end+1} = loaded.data;
    end
    % for each session, compute the expected response timeseries
    % keep only control sessions
    outs = {};
    
    tSeries25 = [];
    for di = 1:length(datas)
        if ~any(datas{di}.design(:,9)==0)
            out = struct;
            
            out.filename = files(di).name;

            clear tSeries25_l tSeries25_r
            % reduce l/r ROIs to just main ROIs
            reduce = [1 2
                      3 4
                      5 6
                      7 8
                      9 10
                      11 12
                      13 14
                      15 16];
            for ri = 1:8
                tSeries25_l(:,ri) = datas{di}.rtSeries25{reduce(ri,1)};
                tSeries25_r(:,ri) = datas{di}.rtSeries25{reduce(ri,2)};
            end
            
            tSeries25 = tSeries25_r-tSeries25_l;

            % pull betas
            for ri = 1:8
                betas(ri,1) = nfit.params.(sprintf('beta_control_%s_cohw',ROIs{ri}));
                betas(ri,2) = nfit.params.(sprintf('beta_control_%s_conw',ROIs{ri}));
            end

            clear tSeries_att

            
            % weight by contrast and coherence
            tSeries_att(:,1) = tSeries25 * betas(:,1);
            tSeries_att(:,2) = tSeries25 * betas(:,2);
            
            % normalize
            for i=1:2
                tSeries_att(:,i) = tSeries_att(:,i) / max(tSeries_att(:,i));
            end
            
            % figure out SV correspondence
            sv_task = zeros(size(tSeries_att,1),2);
            runs = datas{di}.runtrans;
            for ri = 1:size(runs,1)
                rmin = runs(ri,1); rmax = runs(ri,2);
                idxs = logical((datas{di}.design(:,1)>=rmin).*(datas{di}.design(:,1)<=rmax));
                val = mean(datas{di}.design(idxs,9));
                sv_task(rmin:rmax,val) = 1;
            end
            
            out.tSeries = zeros(size(sv_task,1),1);
            for si = 1:size(sv_task,1)
                out.tSeries(si) = tSeries_att(si,find(sv_task(si,:),1));
            end

    %         all = [tS_con tS_coh tSeries25];
    %         
    %         figure; plot(all);

            % the output timeseries will be the result of checking for each
            % dataset 

            outs{end+1} = out;
        end
    end
    aouts{ni} = outs;
end

%% Load relevant files and compute overlay

for ni = 1:length(nSIDs)
    for oi = 1:length(aouts{ni})
        respCorrOverlay(aouts{ni}{oi});
    end
end

%% 
for ni = 1:length(nSIDs)
    for oi = 1:length(aouts{ni})
        view_rca(aouts{ni}{oi});
    end
end
