%% To start off, let's get files

mglSetSID('s300')

global analysis

[analysis.datFolder, analysis.anFolder] = getSubjDataFolder('cohCon');

year = date;
year = year(end-1:end);

files = dir(sprintf('%s/%s*mat',analysis.datFolder,year));

%% Loop over files, fixing them for analysis

fixFiles(files);

%% Loop over files and load them into expHolder

expHolder = loadExp(files);

%% Loop over exp files and export them to csv

for ei = 1:length(expHolder)
%     try
        cohCon2csv(expHolder{ei});
%     catch
%         disp(sprintf('Experiment file %i not generated...',ei));
%     end
end

%% The following functions use the most recent stimulus file

load(fullfile(analysis.datFolder,files(end).name));

%% Generate staircase graphs

% staircaseplots(stimulus);

%% Send main to CSV file
plotting = fade_discFuncs(stimulus,0);
fade_plo2csv(plotting);

%% Send catch to CSV file
cat = fade_catPerf(stimulus);
per2csv(peripheral);