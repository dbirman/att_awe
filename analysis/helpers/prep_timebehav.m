
%% Plot behavioral results (without model fits)

contrasts = [0.325 0.85]-0.25;
coherences = [0.15 0.55];

thresholds = nan(length(sids),2,2,3);
for si = 1:length(sids)
    subj = sids{si};
    files = dir(fullfile(datafolder,'time',subj));
    load(fullfile(datafolder,'time',sprintf('%s',subj),files(end).name));
    staircases = stimulus.staircases.nocatch; 
    % staircases are task | pedestal | length (250/500/1000
    
    for ti = 1:2
        for pi = 1:2
            for li = 1:3
                if staircases{ti,pi,li}(end).trialNum>40
                    out = doStaircase('threshold',staircases{ti,pi,li}(end),'type=weibull');
                else
                    out = doStaircase('threshold',staircases{ti,pi,li},'type=weibull');
                end
                thresholds(si,ti,pi,li) = out.threshold;
            end
        end
    end
end

thresholds(thresholds>1) = NaN;
thresholds(thresholds<0) = NaN;

%% PreP
thresh = squeeze(bootci(1000,@nanmean,thresholds));
thresh_ = squeeze(mean(thresh));