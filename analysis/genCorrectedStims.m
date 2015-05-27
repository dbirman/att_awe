function genCorrectedStims(folder)

%% Get stim files from folder

files = dir(fullfile(folder,'Etc','*.mat'));

%% Mkdir

if ~isdir(fullfile(folder,'Etc','Orig'))
    mkdir(fullfile(folder,'Etc','Orig'));
end

%% Run for each
for i = 1:length(files)
    modifyStim(fullfile(folder,'Etc',files(i).name),fullfile(folder,'Etc','Orig',files(i).name));
end

function modifyStim(file,backup)

%% Load stim

load(file);

copyfile(file,backup);

%% Track file loc

%% Check for vals

if ~exist('stimulus','var') || ~exist('task','var') || ~exist('myscreen','var')
    disp('Load failed for this stimfile, check by hand?');
    keyboard
end

%% Get current stimulus values

randVars = task{1}{1}.randVars;
trials = task{1}{1}.trialnum;
param = task{1}{1}.block.parameter;

baseCoh = randVars.coherence(1:trials);
baseCon = randVars.contrast(1:trials);

lCon = []; lCoh = [];
rCon = []; rCoh = [];

for t = 1:trials
    if param.conSide(t) == 1
        % contrast left
        lCon(t) = baseCon(t) + randVars.conDelta(t);
        rCon(t) = baseCon(t);
    else
        rCon(t) = baseCon(t) + randVars.conDelta(t);
        lCon(t) = baseCon(t);
    end
    if param.cohSide(t) == 1
        lCoh(t) = baseCoh(t) + randVars.cohDelta(t);
        rCoh(t) = baseCoh(t);
    else
        lCoh(t) = baseCoh(t);
        rCoh(t) = baseCoh(t) + randVars.cohDelta(t);
    end
end
lCon = nearestFive(lCon); rCon = nearestFive(rCon);
lCoh = nearestFive(lCoh); rCoh = nearestFive(rCoh);

addCalculatedVar('lCon',lCon,file);
addCalculatedVar('lCoh',lCoh,file);
addCalculatedVar('rCon',rCon,file);
addCalculatedVar('rCoh',rCoh,file);

%% Save

% save(file,'task', 'stimulus', 'myscreen','-V6');

function val = nearestFive(val)

if length(val)>1
    for i = 1:length(val)
        val(i) = nearestFive(val(i));
    end
    return
end

if val > 1 || val < 0
    disp('only works on percentages');
    keyboard
end

val = val*100 - mod(val*100,10) + round(mod(val*100,10)/10)*10;
val = val / 100;
% 
% val = val*100 - mod(val*100,5) + round(mod(val*100,5)/5)*5;
% val = val / 100;