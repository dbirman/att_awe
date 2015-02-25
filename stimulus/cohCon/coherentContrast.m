
% flowAwe
%
%      usage: myscreen=flowAwe()
%         by: daniel birman
%       date: 11/10/14
%    purpose: contrast change detection with cued selective attention.
%
%        use: call flowAwe() to initialize. The first
%             run takes significantly longer due to loading stimuli.
%
%      flags: stimFileNum (-1/#) - Load a specific stimfile from a
%             subject's folder. Defaults to the last stimfile. Warning:
%             Only the first stimfile is saved with the image file data,
%             subsequent stimfiles only contain run data.
%             plots (0/1) - Displays staircase plots (and estimated
%             psychophysic functions)
%             practice (def=0) - Displays on ten trials
% 

function [myscreen] = coherentContrast(varargin)

global stimulus
global fixStimulus
%% Initialize Variables

% add arguments later
stimFileNum = [];
unattended = [];
plots = [];
overrideTask = [];
getArgs(varargin,{'stimFileNum=-1','unattended=0', ...
    'dual=0','plots=1','overrideTask=0'});

stimulus.unattended = unattended;

stimulus.counter = 1; % This keeps track of what "run" we are on.
%% Setup Screen

screen.screenNumber = 2;
myscreen = initScreen(screen);

%% Open Old Stimfile
stimulus.initStair = 1;

if ~isempty(mglGetSID) && isdir(sprintf('~/data/coherentContrast/%s',mglGetSID))
    % Directory exists, check for a stimefile
    files = dir(sprintf('~/data/coherentContrast/%s/1*mat',mglGetSID));

    if length(files) >= 1
        if stimFileNum == -1
            if length(files) > 1
                warning('Multiple stimfiles found, loading last one. If you wanted a different functionality use stimFileNum=#');
            end
            fname = files(end).name;
        else
            fname = files(stimFileNum).name;
        end
        s = load(sprintf('~/data/coherentContrast/%s/%s',mglGetSID,fname));
        stimulus.staircase = s.stimulus.staircase;
        stimulus.stairCatch = s.stimulus.stairCatch;
        stimulus.counter = s.stimulus.counter + 1;

        % load blocks too
        stimulus.runs = s.stimulus.runs;
        stimulus.runs.loaded = 1;

        clear s;
        stimulus.initStair = 0;
        disp(sprintf('(cohCon) Data file: %s loaded, this is run #%i',fname,stimulus.counter));
    end
end

%% Initialize Stimulus

myscreen = initStimulus('stimulus',myscreen);

stimulus.responseKeys = [9 10]; % corresponds to LEFT - RIGHT

%% Colors
stimulus.colors.rmed = 127.75;

% We're going to add an equal number of reserved colors to the top and
% bottom, to try to keep the center of the gamma table stable.
stimulus.colors.reservedBottom = [1 1 1; 0 0 0]; % fixation cross colors
stimulus.colors.reservedTop = [1 0 0; 0 1 0]; % correct/incorrect colors
stimulus.colors.black = 1/255; stimulus.colors.white = 0/255;
stimulus.colors.red = 254/255; stimulus.colors.green = 255/255;
stimulus.colors.nReserved = 2; % this is /2 the true number, because it's duplicated
stimulus.colors.nUnreserved = 256-(2*stimulus.colors.nReserved);

stimulus.colors.mrmax = stimulus.colors.nReserved - 1 + stimulus.colors.nUnreserved;
stimulus.colors.mrmin = stimulus.colors.nReserved;

%% Everything else
stimulus.dots.xcenter = 0;
stimulus.dots.ycenter = 0;
stimulus.dots.dotsize = 3;
stimulus.dots.density = 5;
stimulus.dots.speed = 3;
stimulus.dots.centerOffset = 2;
stimulus.dotsR = stimulus.dots;
stimulus.dotsR.mult = 1;
stimulus.dotsL = stimulus.dots;
stimulus.dotsL.mult = -1;
stimulus = rmfield(stimulus,'dots');

stimulus.pedestals.pedOpts = {'coherence','contrast'};
stimulus.pedestals.coherence = [.1 .3 .5 .7];
stimulus.pedestals.initThresh.coherence = .2;
stimulus.pedestals.contrast = exp(-1.5:(1.25/3):-.25);
stimulus.pedestals.initThresh.contrast = .1;

stimulus.dotsR = initDotsRadial(stimulus.dotsR,myscreen);
stimulus.dotsL = initDotsRadial(stimulus.dotsL,myscreen);

%% Gamma Table Initialization

% get gamma table
if ~isfield(myscreen,'gammaTable')
  stimulus.linearizedGammaTable = mglGetGammaTable;
  disp(sprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'));
  disp(sprintf('(cuecon:initGratings) No gamma table found in myscreen. Contrast displays like this'));
  disp(sprintf('         should be run with a valid calibration made by moncalib for this monitor.'));
  disp(sprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'));
end
stimulus.linearizedGammaTable = myscreen.initScreenGammaTable;
    
%% Character textures
mglTextSet('Helvetica',32,stimulus.colors.black,0,0,0,0,0,0,0);
stimulus.text.mTexK = mglText('M');
stimulus.text.cTexK = mglText('C');
mglTextSet('Helvetica',32,stimulus.colors.white,0,0,0,0,0,0,0);
stimulus.text.mTexW = mglText('M');
stimulus.text.cTexW = mglText('C');

%% MGL Text Parameters
% mglTextSet('Helvetica',32,stimulus.colors.white,0,0,0,0,0,0,0);

%% Setup Task

% This is the contrast change detection task
task{1}{1}.waitForBacktick = 1;
stimulus.seg.ITI = 1; % the ITI is either 20s (first time) or 1s
stimulus.seg.stim = 2;
stimulus.seg.ISI = 3;
stimulus.seg.resp = 4;
task{1}{1}.segmin = [.4 .6 .1 1];
task{1}{1}.segmax = [.8 .6 .4 1];
task{1}{1}.synchToVol = [0 0 0 0];
task{1}{1}.getResponse = [0 0 0 1];
task{1}{1}.parameter.side = [1 2]; % 1 = left, 2 = right, the side will be the one with con/flow + delta (From staircase)
task{1}{1}.parameter.dir = [-1 1];
task{1}{1}.parameter.conPedestal = [1 2 3 4]; % target contrast
task{1}{1}.parameter.cohPedestal = [1 2 3 4]; % target flow coherence
task{1}{1}.parameter.catch = [1 0 0 0 0 0 0 0 0 0]; % 10% chance of being a catch trial
task{1}{1}.random = 1;
task{1}{1}.numTrials = 150;

%% Run variables
task{1}{1}.randVars.calculated.task = nan; % Current task (calc per run)
task{1}{1}.randVars.calculated.deltaPed = nan; % Current task (calc per run)
task{1}{1}.randVars.calculated.coherence = nan;
task{1}{1}.randVars.calculated.contrast = nan;
task{1}{1}.randVars.calculated.trialNum = nan;

%% Tracking

% these are variables that we want to track for later analysis.
task{1}{1}.randVars.calculated.correct = nan;
task{1}{1}.randVars.calculated.trialNum = nan;

stimulus.curTrial = 0;

%% Block setup
if isfield(stimulus,'runs') && isfield(stimulus.runs,'loaded')
    % We already have our blocks
    stimulus.runs = rmfield(stimulus.runs,'loaded'); % remove the load field, otherwise it gets saved across runs
    if stimulus.counter > length(stimulus.runs.taskList)
        stimulus.runs.taskList  = repmat(stimulus.runs.taskList,1,2);
    end
else
    % This is the first run, build up the blocks.
    stimulus.runs = struct;
    stimulus.runs.taskOpts = [1 2];
    stimulus.runs.taskOptsText = {'Motion','Contrast'};
    stimulus.runs.taskList = stimulus.runs.taskOpts(randperm(2));
end
if overrideTask > 0
    stimulus.runs.curTask = overrideTask;
else
    stimulus.runs.curTask = stimulus.runs.taskList(stimulus.counter);
end

%% Unattended Mode
if stimulus.unattended
    fixStimulus.diskSize = 0;
    fixStimulus.stimColor = [.5 .5 .5];
    fixStimulus.responseColor = stimulus.colors.white;
    fixStimulus.interColor = stimulus.colors.black;
    fixStimulus.correctColor = stimulus.colors.green;
    fixStimulus.incorrectColor = stimulus.colors.red;
    [task{2}, myscreen] = fixStairInitTask(myscreen);
end

%% Full Setup
% Initialize task (note phase == 1)
for phaseNum = 1:length(task{1})
    [task{1}{phaseNum}, myscreen] = initTask(task{1}{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@getResponseCallback,@startTrialCallback,[],[]);
end

mglClearScreen(0.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% init staircase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if stimulus.initStair
    % We are starting our staircases from scratch
    disp(sprintf('(cohCon) Initializing staircases'));
    stimulus = initStaircase(stimulus);
else
    disp('(cohCon) Re-using staircase from previous run...');
    % Reset staircase if necessary
    checkStaircaseStop(stimulus);
end

%% EYE CALIB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the eye calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% myscreen = eyeCalibDisp(myscreen);

%% Get Ready...
% clear screen
mglClearScreen(0.5);
mglTextDraw(stimulus.runs.taskOptsText{stimulus.runs.curTask},[0 0]);

% let the user know
disp(sprintf('(cohCon) Starting run number: %i',stimulus.counter));
myscreen.flushMode = 1;

%% Main Task Loop

phaseNum = 1;
% Again, only one phase.
while (phaseNum <= length(task{1})) && ~myscreen.userHitEsc
    % update the task
    [task{1}, myscreen, phaseNum] = updateTask(task{1},myscreen,phaseNum);
  % update the fixation task
    if unattended
        [task{2}, myscreen] = updateTask(task{2},myscreen,1);
    end
    % flip screen
    myscreen = tickScreen(myscreen,task);
end

% task ended
mglClearScreen(0.5);
mglTextDraw('Run complete... please wait.',[0 0]);
mglFlush
myscreen.flushMode = 1;
disp('(cohCon) Run ending...');

if plots
    disp('(cohCon) Displaying plots');
    dispStaircase(stimulus);
end

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);

%%%%%%%%%%%%%%%%%%%%%%%%% EXPERIMENT OVER: HELPER FUNCTIONS FOLLOW %%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Runs at the start of each Trial %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = startTrialCallback(task,myscreen)

global stimulus

stimulus.curTrial = stimulus.curTrial + 1;

%  Set the current task
if stimulus.unattended
    task.thistrial.task = 3;
else
    if task.thistrial.catch
        switchTasks = [2 1];
        task.thistrial.task = switchTasks(stimulus.runs.curTask);
        % edit seglen
        task.thistrial.seglen(stimulus.seg.ISI) = .5;
        task.thistrial.seglen(stimulus.seg.resp) = 2;
        disp('(cohCon) Catch trial.');
    else
        task.thistrial.task = stimulus.runs.curTask;
    end
end

% Set the missing thistrial vars
task.thistrial.coherence = stimulus.pedestals.coherence(task.thistrial.cohPedestal);
task.thistrial.contrast = stimulus.pedestals.contrast(task.thistrial.conPedestal);
task.thistrial.trialNum = stimulus.curTrial;
[task.thistrial.deltaPed, stimulus] = getDeltaPed(task,stimulus,task.thistrial.task);

% Assign the deltaPed to the correct locations
if task.thistrial.task==1
    % speed
    stimulus.live.cohDelta = task.thistrial.deltaPed;
    stimulus.live.conDelta = 0;
    disp(sprintf('(cohCon) Trial %i starting. Coherence: %.02f + %.02f Contrast %.02f',task.thistrial.trialNum,task.thistrial.coherence,stimulus.live.cohDelta,task.thistrial.contrast));
elseif task.thistrial.task==2
    % contrast
    stimulus.live.cohDelta = 0;
    stimulus.live.conDelta = task.thistrial.deltaPed;
    if (task.thistrial.contrast + stimulus.live.conDelta) > 1
        stimulus.live.conDelta = 1 - task.thistrial.contrast;
    end
    disp(sprintf('(cohCon) Trial %i starting. Coherence: %.02f Contrast %.02f + %.02f',task.thistrial.trialNum,task.thistrial.coherence,task.thistrial.contrast,stimulus.live.conDelta));
else
    % unattended
    stimulus.live.cohDelta = 0;
    stimulus.live.conDelta = 0;
    disp(sprintf('(cohCon) Trial %i starting. Coherence: %.02f Contrast %.02f',task.thistrial.trialNum,task.thistrial.coherence,task.thistrial.contrast));
end

% set the gammaTable for this trial
if ~stimulus.unattended
    setGammaTable_flowMax(task.thistrial.contrast + stimulus.live.conDelta);
else
    setGammaTable_flowMax(1);
end

% set directions
stimulus.dotsL.dir = task.thistrial.dir;
stimulus.dotsR.dir = task.thistrial.dir;

function ped = curPedValue(task)
if task.thistrial.task==1
    ped = task.thistrial.cohPedestal;
else
    ped = task.thistrial.conPedestal;
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Runs at the start of each Segment %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = startSegmentCallback(task, myscreen)

myscreen.flushMode = 0;

global stimulus

switch task.thistrial.thisseg
    case stimulus.seg.ITI
        stimulus.live.dots = 0;
        stimulus.live.fixColor = stimulus.colors.black;
        stimulus.live.catchFix = 0;
    case stimulus.seg.stim
        stimulus.live.dots = 1;
        stimulus.live.fixColor = stimulus.colors.black;
        stimulus.live.catchFix = 0;
    case stimulus.seg.ISI
        stimulus.live.dots = 0;
        stimulus.live.fixColor = stimulus.colors.black;
        stimulus.live.catchFix = 1;
    case stimulus.seg.resp
        stimulus.live.dots = 0;
        stimulus.live.fixColor = stimulus.colors.white;
        stimulus.live.catchFix = 1;
end
if task.thistrial.thisseg == stimulus.seg.stim
    stimulus.live.dots = 1;
else
    stimulus.live.dots = 0;
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Refreshes the Screen %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = screenUpdateCallback(task, myscreen)
global stimulus

mglClearScreen(0.5);

if stimulus.live.dots==1, stimulus = upDots(task,stimulus,myscreen); end
if ~stimulus.unattended, upFix(task,stimulus); end

%%
function upFix(task,stimulus)

if ~task.thistrial.catch || stimulus.live.catchFix == 0
    mglFixationCross(1,1,stimulus.live.fixColor);
else
    if task.thistrial.task==1
        if stimulus.live.fixColor==stimulus.colors.white
            mglBltTexture(stimulus.text.mTexW,[0 0]);
        else
            mglBltTexture(stimulus.text.mTexK,[0 0]);
        end
    else
        if stimulus.live.fixColor==stimulus.colors.white
            mglBltTexture(stimulus.text.cTexW,[0 0]);
        else
            mglBltTexture(stimulus.text.cTexK,[0 0]);
        end
    end
end

function stimulus = upDots(task,stimulus,myscreen)

% update the dots

if task.thistrial.side == 1
    stimulus.dotsL = updateDotsRadial(stimulus.dotsL,task.thistrial.coherence+stimulus.live.cohDelta,myscreen);
    stimulus.dotsR = updateDotsRadial(stimulus.dotsR,task.thistrial.coherence,myscreen);
else
    stimulus.dotsL = updateDotsRadial(stimulus.dotsL,task.thistrial.coherence,myscreen);
    stimulus.dotsR = updateDotsRadial(stimulus.dotsR,task.thistrial.coherence+stimulus.live.cohDelta,myscreen);
end

if task.thistrial.side == 1
    lConDelta = stimulus.live.conDelta;
    rConDelta = 0;
else
    lConDelta = 0;
    rConDelta = stimulus.live.conDelta;
end

% Correct values for gamma table adjustments
correctCon = task.thistrial.contrast / stimulus.curMaxContrast;
rConDelta = rConDelta / stimulus.curMaxContrast;
lConDelta = lConDelta / stimulus.curMaxContrast;

% Correct values for size of gamma table

% dotsR
% update +contrast

mglPoints2(stimulus.dotsR.xdisp(stimulus.dotsR.con==1),stimulus.dotsR.ydisp(stimulus.dotsR.con==1),...
    stimulus.dotsR.dotsize,[.5 .5 .5] - adjustConToTable(correctCon + rConDelta,stimulus)/2);
% update - contrast
mglPoints2(stimulus.dotsR.xdisp(stimulus.dotsR.con==2),stimulus.dotsR.ydisp(stimulus.dotsR.con==2),...
    stimulus.dotsR.dotsize,[.5 .5 .5] + adjustConToTable(correctCon + rConDelta,stimulus)/2);
% dotsL
% update +contrast
mglPoints2(stimulus.dotsL.xdisp(stimulus.dotsL.con==1),stimulus.dotsL.ydisp(stimulus.dotsL.con==1),...
    stimulus.dotsL.dotsize,[.5 .5 .5] - adjustConToTable(correctCon + lConDelta,stimulus)/2);
% update - contrast
mglPoints2(stimulus.dotsL.xdisp(stimulus.dotsL.con==2),stimulus.dotsL.ydisp(stimulus.dotsL.con==2),...
    stimulus.dotsL.dotsize,[.5 .5 .5] + adjustConToTable(correctCon + lConDelta,stimulus)/2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Adjust contrast to the gamma table %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function conValue = adjustConToTable(conValue,stimulus)
conValue = conValue * stimulus.colors.nUnreserved / 256;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Called When a Response Occurs %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function [task, myscreen] = getResponseCallback(task, myscreen)

global stimulus

responseText = {'Incorrect','Correct'};
fixColors = {stimulus.colors.red,stimulus.colors.green};

if any(task.thistrial.whichButton == stimulus.responseKeys)
    if task.thistrial.gotResponse == 0
        task.thistrial.correct = task.thistrial.whichButton == stimulus.responseKeys(task.thistrial.side);
        % Store whether this was correct
        stimulus.live.fixColor = fixColors{task.thistrial.correct+1};
        disp(sprintf('(cohCon) Response %s',responseText{task.thistrial.correct+1}));
        if ~task.thistrial.catch
            stimulus.staircase{task.thistrial.task,curPedValue(task)} = ...
                doStaircase('update',stimulus.staircase{task.thistrial.task,curPedValue(task)},task.thistrial.correct);
        else
            stimulus.live.fixColor = stimulus.colors.black; % we never show information about catch trials
            stimulus.live.catchFix = 0;
            stimulus.stairCatch{task.thistrial.task} = ...
                doStaircase('update',stimulus.stairCatch{task.thistrial.task},task.thistrial.correct);
        end
    else
        disp(sprintf('(cohCon) Subject responded multiple times: %i',task.thistrial.gotResponse+1));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                              HELPER FUNCTIONS                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% getDeltaPed

function [deltaPed, stimulus] = getDeltaPed(task,stimulus,taskNum)
if task.thistrial.catch
    [deltaPed, stimulus.stairCatch{taskNum}] = doStaircase('testValue',stimulus.stairCatch{taskNum});
else
    [deltaPed, stimulus.staircase{taskNum,curPedValue(task)}] = doStaircase('testValue',stimulus.staircase{taskNum,curPedValue(task)});
end


%%%%%%%%%%%%%%%%%%%%%%%%
%    initStaircase     %
%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initStaircase(stimulus)

stimulus.stairCatch = cell(1,2);
stimulus.staircase = cell(2,length(stimulus.pedestals.contrast));
stimulus.stairCatch{1} = doStaircase('init','fixed',...
    'fixedVals',[.1 .14 .18 .22 .26]);
stimulus.stairCatch{2} = doStaircase('init','fixed',...
    'fixedVals',[.05 .075 .1 .125 .15]);
for task = 1:2
    stimulus.staircase{task,1} = doStaircase('init','upDown',...
        'initialThreshold',stimulus.pedestals.initThresh.(stimulus.pedestals.pedOpts{task}),...
        'initialStepsize',stimulus.pedestals.initThresh.(stimulus.pedestals.pedOpts{task})/3,...
        'minThreshold=0.001','maxThreshold=0.5','stepRule','pest', ...
        'nTrials=50','maxStepsize=.2','minStepsize=.001');
end

for i = 2:length(stimulus.pedestals.coherence)
    stimulus.staircase{1,i} = stimulus.staircase{1,1};
    stimulus.staircase{2,i} = stimulus.staircase{2,1};
end

%%%%%%%%%%%%%%%%%%%%%%%
%    dispStaircase    %
%%%%%%%%%%%%%%%%%%%%%%%
function dispStaircase(stimulus)

try
    taskOpts = {'coherence','contrast'};
    
    plotting = zeros(2,4);
    
    drawing = {'-r' '-g' '-b' '-y'
                '--r' '--g' '--b' '--y'};
    for task = 1:2
        figure % this is the 'staircase' figure
        title(sprintf('%s, Staircase plot (R->G->B->Y high)',taskOpts{task}));
        hold on
        for ped = 1:4
            try
                testV = [];
                for i = 1:length(stimulus.staircase{task,ped})
                    testV = [testV stimulus.staircase{task,ped}(i).testValues];
                end
                plot(testV,drawing{task,ped});
            catch
            end
            try
                out = doStaircase('threshold',stimulus.staircase{task,ped},'type','weibull'); % noise, 1 cue, lowest
                plotting(task,ped) = out.threshold;
            catch
                plotting(task,ped) = -1;
            end
        end
    end
    hold off
    figure
    hold on
    title(sprintf('%s, R->G->B High',taskOpts{task}));
    plot(stimulus.pedestals.(taskOpts{task})(1:4),plotting(1,:),'-r');
    plot(stimulus.pedestals.(taskOpts{task})(1:4),plotting(2,:),'--r');
    legend(taskOpts);
    axis([stimulus.pedestals.(taskOpts{task})(1) stimulus.pedestals.(taskOpts{task})(4) 0 .3]);
    hold off

catch
    disp('(cohCon) Figures were not generated successfully.');
end

%% checkStaircaseStop
function checkStaircaseStop(stimulus)
% Check both staircases
for task = 1:2
    for ped = 1:4
        s = stimulus.staircase{task,ped};
        if doStaircase('stop',s)
            % this is a bit of a pain... you can't pass an initialThreshold
            % argument do doStaircase('init',s, ...), it ignores everything and
            % resets using the calculated threshold. Because you can't override it
            [args, vals, ~] = getArgs(s(1).initArgs);
            threshPos = -1;
            stepPos = -1;
            for i = 1:length(args)
                switch args{i}
                    case 'initialThreshold'
                        threshPos = i;
                    case 'initialStepsize'
                        stepPos = i;
                end
            end
            out = doStaircase('threshold',s);
            in = input(sprintf('Resetting Staircase... Estimate is: %1.2f. Reset ([Y]/[C]ustom/[O]riginal): ',out.threshold),'s');
            switch in
                case 'Y'
                    vals{threshPos} = out.threshold;
                    vals{stepPos} = out.threshold / 3;
                case 'C'
                    disp('Original values:');
                    disp(sprintf('%s: %0.2f',args{threshPos},num2str(vals{threshPos})));
                    val = str2double(input('New threshold value: ','s'));
                    vals{threshPos} = val;
                    vals{stepPos} = val / 3;
                case 'O'
            end
<<<<<<< HEAD
            stimulus.staircase{task,ped}(end+1) = doStaircase('init',s,'initialThreshold',vals{threshPos},'initialStepsize',vals{stepPos});
=======
            stimulus.staircase{task,ped} = doStaircase('init',s,args{1},vals{2},args{2},vals{2},{args{3},vals{3},...
                args{4},vals{4},args{5},vals{5},args{6},vals{6},args{7},vals{7},args{8},vals{8});
%             stimulus.staircase{task,ped} = doStaircase('init',s,'initialThreshold',vals{threshPos},'initialStepsize',vals{stepPos});
>>>>>>> origin/master
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create dots for optic flow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = initDotsRadial(dots,~)

% maximum depth of points
dots.minX = 3;
dots.maxX = 13;
dots.minY = -5;
dots.maxY = 5;

% make a some points
dots.n = 500*dots.density;
% make sure it's an even number
dots.n = dots.n + mod(dots.n,2);

% set half to white and half to black
dots.con = repmat([1 2],1,dots.n/2);

dots.x = rand(1,dots.n)*(dots.maxX-dots.minX)+dots.minX;
dots.y = rand(1,dots.n)*abs(dots.maxY-dots.minY)+dots.minY;

dots.xdisp = dots.mult*dots.x;
dots.ydisp = dots.y;

% set incoherent dots to 0
dots.coherency = 1;
dots.incoherent = rand(1,dots.n) > dots.coherency;
dots.incoherentn = sum(dots.incoherent);
dots.coherent = ~dots.incoherent;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step dots for Radial
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = updateDotsRadial(dots,coherence,myscreen)
% stuff to compute median speed
dots.oldx = dots.x;
dots.oldy = dots.y;

% get the coherent and incoherent dots
% if (dots.coherency ~= dots.coherence)
  dots.incoherent = rand(1,dots.n) > coherence;
  dots.incoherentn = sum(dots.incoherent);
  dots.coherent = ~dots.incoherent;
  dots.coherency = coherence;
% end
freq_factor = dots.speed/myscreen.framesPerSecond;

% move coherent dots
dots.x(dots.coherent) = dots.x(dots.coherent) + dots.dir*freq_factor;

% these are for flipping into the other quadrants
xmat = [1 1 -1 -1];
ymat = [1 -1 1 -1];
xmat = repmat(xmat,1,dots.incoherentn+4-mod(dots.incoherentn,4));
ymat = repmat(ymat,1,dots.incoherentn+4-mod(dots.incoherentn,4));
xmat = xmat(randperm(dots.incoherentn));
ymat = ymat(randperm(dots.incoherentn));

% move incoherent dots
% get random vectors
dots.rX = rand(1,dots.incoherentn); % get a random number -1 -> 1 for the x offset
dots.rY = sqrt(1-dots.rX.^2); % solve for y
dots.rX = dots.rX .* xmat;
dots.rY = dots.rY .* ymat;
dots.rX = dots.rX * freq_factor; % rescale to match the velocity
dots.rY = dots.rY * freq_factor;
dots.x(dots.incoherent) = dots.x(dots.incoherent) + dots.rX;
dots.y(dots.incoherent) = dots.y(dots.incoherent) + dots.rY;

offscreen = dots.x > dots.maxX;
dots.x(offscreen) = dots.x(offscreen) - abs(dots.maxX - dots.minX);
offscreen = dots.x < dots.minX;
dots.x(offscreen) = dots.x(offscreen) + abs(dots.maxX - dots.minX);

offscreen = dots.y > dots.maxY;
dots.y(offscreen) = dots.y(offscreen) - abs(dots.maxY - dots.minY);
offscreen = dots.y < dots.minY;
dots.y(offscreen) = dots.y(offscreen) + abs(dots.maxY - dots.minY);


dots.xdisp = dots.mult*dots.x;
dots.ydisp = dots.y;