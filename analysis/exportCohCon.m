%% To start off, let's get files

addpath(genpath('~/proj/att_awe/analysis/'));
doEye = true;
scan = false;

global analysis

[analysis.datFolder, analysis.anFolder] = getSubjDataFolder('cohcon',scan);

year = date;
year = year(end-1:end);

files = dir(sprintf('%s/%s*mat',analysis.datFolder,year));

%% Loop over files, fixing them for analysis

fixFiles(files);

%% Loop over files and load them into expHolder

expHolder = loadExp(files);

%% Loop over exp files and export them to csv

for ei = 1:length(expHolder)
    %     try                        %EYE SKIP
    cohCon2csv(expHolder{ei},~doEye);
    %     catch
    %         disp(sprintf('Experiment file %i not generated...',ei));
    %     end
end
if doEye
    return
end

%% The following functions use the most recent stimulus file

load(fullfile(analysis.datFolder,files(end).name));

%% If Scan: Concatenate across stimulus runs
if ~isempty(strfind(analysis.anFolder,'scan'))
    taskCued = {};
    taskMiscued = {};
    control = {};
    for i = 1:length(expHolder)
        expi = expHolder{i}{1};
        for j = 1:2
            taskCued{i,j} = expi.stimulus.staircase{j};
            taskMiscued{i,j} = expi.stimulus.stairCatch{j};
            for k = 1:size(expi.stimulus.nocatchs.staircase,2)
                
%                 try
                control{i,j,k} = expi.stimulus.nocatchs.staircase{j,k};
%                 catch
%                 end
            end
        end
    end
end

%% Send main to CSV file
[plotting, mfits] = cohCon_discFuncs(stimulus,0);
cohCon_plo2csv(plotting);

%% Send catch to CSV file
[cat, cfits] = cohCon_catPerf(stimulus,0);
per2csv(cat);

%% Send nocatch to CSV file
[nocat, nfits] = cohCon_nocatPerf(stimulus,0);
nocat2csv(nocat);

%% Send fits to CSV file
cohCon_fits2csv(nfits,mfits,cfits);

% %%
% x = 0:.01:3
% figure, hold on
% color = {{'k','g','k','k'},{'r','b','r','r'}};
% for t = 1:2
%     dat = mfits{t,1};
%     dat(:,1) = dat(:,1) ./ mean(nfits{t,1}(:,1),1);
%     y = weibull(x,mean(dat,1));
%     plot(x,y,color{t}{1});
% end
% for t = 1:2
%     dat = cfits{t,1};
%     dat(:,1) = dat(:,1) ./ mean(nfits{t,1}(:,1),1);
%     y = weibull(x,mean(dat,1));
%     plot(x,y,color{t}{2});
% end
% axis([0 5 .5 1])

%% normalize
