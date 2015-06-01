function genCorrectedStims(folder)

%% Get stim files from folder

files = dir(fullfile(folder,'Etc','*.mat'));

%% Mkdir

if ~isdir(fullfile(folder,'Etc','Orig'))
    mkdir(fullfile(folder,'Etc','Orig'));
end

%% Run for each
for i = 1:length(files)
    try
    modifyStim(fullfile(folder,'Etc',files(i).name),fullfile(folder,'Etc','Orig',files(i).name));
    catch
    end
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

nTask = [];

lCon = []; lCoh = [];
rCon = []; rCoh = [];

if any(param.catch(1:trials)>0)
    nTask = randVars.task(1:trials) + param.catch(1:trials) .* 2;
else
    nTask = randVars.task(1:trials);
end

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

if ~isfield(randVars,'lCon'), addCalculatedVar('lCon',lCon,file,'backup=0'); end
if ~isfield(randVars,'lCoh'), addCalculatedVar('lCoh',lCoh,file,'backup=0'); end
if ~isfield(randVars,'rCon'), addCalculatedVar('rCon',rCon,file,'backup=0'); end
if ~isfield(randVars,'rCoh'), addCalculatedVar('rCoh',rCoh,file,'backup=0'); end
if ~isfield(randVars,'nTask'), addCalculatedVar('nTask',nTask,file,'backup=0'); end

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