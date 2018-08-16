%% Permutation test
% Figure out the range of expected R^2 on random data (take the individual
% subject and randomly permute the responses, column 8)
% reps = 500; % expect to take about 30 minutes per repetition
% 
reps = 500;
permlike = zeros(reps,length(nSIDs),2,2);

tic;
for rep = 1:reps/50
    % do 50 repeats at a time
    localr2 = zeros(50,length(nSIDs),2,2);
    disppercent(-1/50);
    for irep = 1:50
        
        parr2 = localr2(irep,:,:,:);
        for ai = 1:length(nSIDs)
            adata = loadadata(sprintf('s%03.0f',nSIDs(ai)));

            pdata = adata;
            pdata(:,8) = pdata(randperm(size(pdata,1)),8);
            ir2 = parr2(ai,:,:);
            for mi = 1:length(mopts)
                for ni = 1:length(bmodels)
                    info = struct;
                    info.sigma = 1;
                    info.model = bmodels{ni};
                    info.rois = [1 8];
                    info.lapse = lapses(ai);
                    info.respcon = squeeze(respcon(ai,:,:));
                    info.respcoh = squeeze(respcoh(ai,:,:));
                    %                     if clapse==0
                    %                         clapse = min(lapses(lapses>0));
                    %                     end
                    fit = fitCCBehavControlModel_fmri(pdata,info,1);
                    ir2(mi,ni) = fit.likelihood;
                end
                close all
            end
            parr2(ai,:,:) = ir2;
        end
        localr2(irep,:,:,:) = parr2;
        disppercent(irep/50);
    end
    disppercent(inf);
    
    idxs = ((rep-1)*50+1):rep*50;
    permlike(idxs,:,:,:) = localr2;
    save(fullfile(datafolder,'permutation.mat','permlike'));
end
t = toc;
disp(sprintf('Took %01.2f seconds',t));


%% Computation permutation R^2 range

quantile(permr2(:),[0.05 0.95 0.99])
hist(permr2(:));


significance = [.05 .01 .001];
pval = [0.0201601085600963 0.0278872183918559 0.0326540667595963];

%% Compute all possible difference scores
pr2 = permr2(:);
diffs = zeros(length(pr2),length(pr2)-1);

for i = 1:length(pr2)
    diffs(i,:) = pr2(setdiff(1:length(pr2),i))-pr2(i);
end
quantile(diffs(:),[0.01 0.05 0.95 0.99])
hist(diffs(:));

%   Columns 1 through 3
% 
%        -0.0276701816421294       -0.0203526997937312        0.0203526997937312
% 
%   Column 4
% 
%         0.0276701816421294