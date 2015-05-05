%% To start off, let's get files

mglSetSID('s025')
addpath(genpath('~/proj/att_awe/analysis/'));

global analysis

[analysis.datFolder, analysis.anFolder] = getSubjDataFolder('cohcon',false);

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
        cohCon2csv(expHolder{ei},true);
%     catch
%         disp(sprintf('Experiment file %i not generated...',ei));
%     end
end

%% The following functions use the most recent stimulus file

load(fullfile(analysis.datFolder,files(end).name));

%% Generate staircase graphs

% staircaseplots(stimulus);

%% Send main to CSV file
plotting = cohCon_discFuncs(stimulus,0);
cohCon_plo2csv(plotting);

%% Send catch to CSV file
cat = cohCon_catPerf(stimulus,0);
per2csv(cat);

%% Send nocatch to CSV file
nocat = cohCon_nocatPerf(stimulus,0);
nocat2csv(nocat);