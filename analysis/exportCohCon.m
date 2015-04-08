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
    try
        cohCon2csv(expHolder{ei});
    catch
        disp(sprintf('Experiment file %i not generated...',ei));
    end
end

%% The following functions use the most recent stimulus file

load(fullfile(analysis.datFolder,files(end).name));

%% Generate staircase graphs

% staircaseplots(stimulus);

%% Send main to CSV file
plotting = cohCon_discFuncs(stimulus,0);
cohCon_plo2csv(plotting);

%% Send catch to CSV file
% temporary code cause I messed up
s1 = stimulus.stairCatch{1};
s2 = stimulus.stairCatch{2};
if size(stimulus.stairCatch,2)~=5
    stimulus.stairCatch = cell(2,5);
    for i = 1:5
        stimulus.stairCatch{1,i} = s1;
        stimulus.stairCatch{2,i} = s2;
    end
end
cat = cohCon_catPerf(stimulus,0);
per2csv(cat);