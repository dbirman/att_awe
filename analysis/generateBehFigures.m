%% To start off, let's get files

mglSetSID('s300')

global analysis

[analysis.datFolder, analysis.anFolder] = getSubjDataFolder;

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
        exp2csv(expHolder{ei});
    catch
        disp(sprintf('Experiment file %i not generated...',ei));
    end
end

%% The following functions use the most recent stimulus file

load(fullfile(analysis.datFolder,files(end).name));

%% Generate staircase graphs

staircaseplots(stimulus);

%% Send main to CSV file
plotting = gen_discFuncs(stimulus,0);
plo2csv(plotting);

%% Send peripheral to CSV file
gen_perPerf(stimulus);
per2csv(peripheral);

%% Eye Movement Analysis

% First find the expHolder files with eye data
eyeData = [];
for e = 1:length(expHolder)
    if length(expHolder{e}) > 3
        eyeData = [eyeData e];
    end
end
eyeHolder = expHolder(eyeData);

%% Analyze

for e = 1:length(eyeHolder)
    eye = eyeHolder{e};
    dispEyeTracesMain(eye);
end