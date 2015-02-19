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

function [myscreen] = flowAwe(varargin)

global stimulus
global fixStimulus
%% Initialize Variables

% add arguments later
stimFileNum = [];
unattended = [];
plots = [];
overrideTask = [];
getArgs(varargin,{'stimFileNum=-1','unattended=0', ...
    'dual=0','plots=1','practice=0','overrideTask=0'});

stimulus.unattended = unattended;

stimulus.counter = 1; % This keeps track of what "run" we are on.
%% Setup Screen

screen.screenNumber = 2;
myscreen = initScreen(screen);

%% Open Old Stimfile
stimulus.initStair = 1;

if ~isempty(mglGetSID) && isdir(sprintf('~/data/flowAwe/%s',mglGetSID))
    % Directory exists, check for a stimefile
    files = dir(sprintf('~/data/flowAwe/%s/1*mat',mglGetSID));

    if length(files) >= 1
        if stimFileNum == -1
            if length(files) > 1
                warning('Multiple stimfiles found, loading last one. If you wanted a different functionality use stimFileNum=#');
            end
            fname = files(end).name;
        else
            fname = files(stimFileNum).name;
        end
        s = load(sprintf('~/data/flowAwe/%s/%s',mglGetSID,fname));
        stimulus.staircase = s.stimulus.staircase;
        stimulus.stairCatch = s.stimulus.stairCatch;
        stimulus.counter = s.stimulus.counter + 1;

        % load blocks too
        stimulus.runs = s.stimulus.runs;
        stimulus.runs.loaded = 1;

        clear s;
        stimulus.initStair = 0;
        disp(sprintf('(flowAwe) Data file: %s loaded, this is run #%i',fname,stimulus.counter));
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
stimulus.dots.dotsize = 4;
stimulus.dots.density = 10;
stimulus.dots.coherence = 1;
stimulus.dots.speed = 10;
stimulus.dots.speedmult = 1;
stimulus.dots.centerOffset = 2;
stimulus.dots.T = [0 0 stimulus.dots.speed/myscreen.framesPerSecond];
stimulus.dotsR = stimulus.dots;
stimulus.dotsR.mult = 1;
stimulus.dotsL = stimulus.dots;
stimulus.dotsL.mult = -1;
stimulus = rmfield(stimulus,'dots');

stimulus.pedestals.pedOpts = {'flow','contrast'};
stimulus.pedestals.flow = [.1 .3 .5 .7];
stimulus.pedestals.initThresh.flow = .2;
stimulus.pedestals.contrast = exp(-1.5:(1.25/3):-.25);
stimulus.pedestals.initThresh.contrast = .05;

stimulus.dotsR = initDotsOpticflow(stimulus.dotsR,myscreen);
stimulus.dotsL = initDotsOpticflow(stimulus.dotsL,myscreen);

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
    
%% Additional Dots Params

% set color
stimulus.dotsR.color = ones(stimulus.dotsR.n,1);
stimulus.dotsL.color = ones(stimulus.dotsL.n,1);

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
task{1}{1}.segmin = [.8 .8 .1 1.1];
task{1}{1}.segmax = [.8 .8 .5 1.1];
task{1}{1}.synchToVol = [0 0 0 0];
task{1}{1}.getResponse = [0 0 0 1];
task{1}{1}.parameter.side = [1 2]; % 1 = left, 2 = right, the side will be the one with con/flow + delta (From staircase)
task{1}{1}.parameter.dir = [-1 1];
task{1}{1}.parameter.conPedestal = [1 2 3 4]; % target contrast
task{1}{1}.parameter.floPedestal = [1 2 3 4]; % target flow coherence
task{1}{1}.parameter.catch = [1 0 0 0 0 0 0 0 0 0]; % 10% chance of being a catch trial
task{1}{1}.random = 1;
task{1}{1}.numTrials = 125;

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
    stimulus.runs.taskOptsText = {'Flow','Contrast'};
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
    disp(sprintf('(flowAwe) Initializing staircases'));
    stimulus = initStaircase(stimulus);
else
    disp('(flowAwe) Re-using staircase from previous run...');
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
disp(sprintf('(flowAwe) Starting run number: %i',stimulus.counter));
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
disp('(flowAwe) Run ending...');

if plots
    disp('(flowAwe) Displaying plots');
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
        disp('(flowAwe) Catch trial.');
    else
        task.thistrial.task = stimulus.runs.curTask;
    end
end

% Set the missing thistrial vars
task.thistrial.coherence = stimulus.pedestals.flow(task.thistrial.floPedestal);
task.thistrial.contrast = stimulus.pedestals.contrast(task.thistrial.conPedestal);
task.thistrial.trialNum = stimulus.curTrial;
[task.thistrial.deltaPed, stimulus] = getDeltaPed(task,stimulus,task.thistrial.task,curPedValue(task));

% Assign the deltaPed to the correct locations
if task.thistrial.task==1
    % coherence
    stimulus.live.cohDelta = task.thistrial.deltaPed;
    if (task.thistrial.coherence + stimulus.live.cohDelta) > .95
        stimulus.live.cohDelta = .95 - task.thistrial.coherence;
    end
    stimulus.live.conDelta = 0;
    disp(sprintf('(flowAwe) Trial %i starting. Coherence: %.02f + %.02f Contrast %.02f',task.thistrial.trialNum,task.thistrial.coherence,stimulus.live.cohDelta,task.thistrial.contrast));
elseif task.thistrial.task==2
    % contrast
    stimulus.live.cohDelta = 0;
    stimulus.live.conDelta = task.thistrial.deltaPed;
    if (task.thistrial.contrast + stimulus.live.conDelta) > 1
        stimulus.live.conDelta = 1 - task.thistrial.contrast;
    end
    disp(sprintf('(flowAwe) Trial %i starting. Coherence: %.02f Contrast %.02f + %.02f',task.thistrial.trialNum,task.thistrial.coherence,task.thistrial.contrast,stimulus.live.conDelta));
else
    % unattended
    stimulus.live.cohDelta = 0;
    stimulus.live.conDelta = 0;
    disp(sprintf('(flowAwe) Trial %i starting. Coherence: %.02f Contrast %.02f',task.thistrial.trialNum,task.thistrial.coherence,task.thistrial.contrast));
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

% randomize speeds
stimulus.dotsL.speedmult = rand*1.5 + 0.5;
stimulus.dotsR.speedmult = rand*1.5 + 0.5;

function ped = curPedValue(task)
if task.thistrial.task==1
    ped = task.thistrial.floPedestal;
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
    stimulus.dotsL = updateDotsOpticFlow(stimulus.dotsL,task.thistrial.coherence+stimulus.live.cohDelta,myscreen);
    stimulus.dotsR = updateDotsOpticFlow(stimulus.dotsR,task.thistrial.coherence,myscreen);
else
    stimulus.dotsL = updateDotsOpticFlow(stimulus.dotsL,task.thistrial.coherence,myscreen);
    stimulus.dotsR = updateDotsOpticFlow(stimulus.dotsR,task.thistrial.coherence+stimulus.live.cohDelta,myscreen);
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

mglPoints2(stimulus.dotsR.x(stimulus.dotsR.con==1),stimulus.dotsR.y(stimulus.dotsR.con==1),...
    stimulus.dotsR.dotsize,[.5 .5 .5] - adjustConToTable(correctCon + rConDelta,stimulus)/2);
% update - contrast
mglPoints2(stimulus.dotsR.x(stimulus.dotsR.con==2),stimulus.dotsR.y(stimulus.dotsR.con==2),...
    stimulus.dotsR.dotsize,[.5 .5 .5] + adjustConToTable(correctCon + rConDelta,stimulus)/2);
% dotsL
% update +contrast
mglPoints2(stimulus.dotsL.x(stimulus.dotsL.con==1),stimulus.dotsL.y(stimulus.dotsL.con==1),...
    stimulus.dotsL.dotsize,[.5 .5 .5] - adjustConToTable(correctCon + lConDelta,stimulus)/2);
% update - contrast
mglPoints2(stimulus.dotsL.x(stimulus.dotsL.con==2),stimulus.dotsL.y(stimulus.dotsL.con==2),...
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
        disp(sprintf('(flowAwe) Response %s',responseText{task.thistrial.correct+1}));
        if ~task.thistrial.catch
            stimulus.staircase{task.thistrial.task,task.thistrial.floPedestal,task.thistrial.conPedestal} = ...
                doStaircase('update',stimulus.staircase{task.thistrial.task,task.thistrial.floPedestal,task.thistrial.conPedestal},task.thistrial.correct);
        else
            stimulus.live.fixColor = stimulus.colors.black; % we never show information about catch trials
            stimulus.live.catchFix = 0;
            stimulus.stairCatch{task.thistrial.task} = ...
                doStaircase('update',stimulus.stairCatch{task.thistrial.task},task.thistrial.correct);
        end
    else
        disp(sprintf('(flowAwe) Subject responded multiple times: %i',task.thistrial.gotResponse+1));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                              HELPER FUNCTIONS                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% getDeltaPed

function [deltaPed, stimulus] = getDeltaPed(task,stimulus,taskNum,pedNum)
if task.thistrial.catch
    [deltaPed, stimulus.stairCatch{taskNum}] = doStaircase('testValue',stimulus.stairCatch{taskNum});
else
    [deltaPed, stimulus.staircase{taskNum,task.thistrial.floPedestal,task.thistrial.conPedestal}] = doStaircase('testValue',stimulus.staircase{taskNum,task.thistrial.floPedestal,task.thistrial.conPedestal});
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
    stimulus.staircase{task,1,1} = doStaircase('init','upDown',...
        'initialThreshold',stimulus.pedestals.initThresh.(stimulus.pedestals.pedOpts{task}),...
        'initialStepsize',stimulus.pedestals.initThresh.(stimulus.pedestals.pedOpts{task})/3,...
        'minThreshold=0.001','maxThreshold=0.5','stepRule','pest', ...
        'nTrials=50','maxStepsize=.2','minStepsize=.001');
end

for i = 1:length(stimulus.pedestals.flow)
    for j = 1:length(stimulus.pedestals.contrast)
        stimulus.staircase{1,i,j} = stimulus.staircase{1,1,1};
        stimulus.staircase{2,i,j} = stimulus.staircase{2,1,1};
    end
end

%%%%%%%%%%%%%%%%%%%%%%%
%    dispStaircase    %
%%%%%%%%%%%%%%%%%%%%%%%
function dispStaircase(stimulus)

try
    taskOpts = {'flow','contrast'};
    
    drawing = {'-r' '-g' '-b' '-y'
                '--r' '--g' '--b' '--y'};
    for task = 1:2
        figure % this is the 'staircase' figure
        title(sprintf('%s, Staircase plot (R->G->B->Y high)',taskOpts{task}));
        hold on
        for fped = 1:4
            for cped = 1:4
                try
                    testV = [];
                    for i = 1:length(stimulus.staircase{task,fped,cped})
                        testV = [testV stimulus.staircase{task,fped,cped}(i).testValues];
                    end
                    plot(testV,drawing{task,ped});
                catch
                end
                try
                    out = doStaircase('threshold',stimulus.staircase{task,fped,cped},'type','weibull'); % noise, 1 cue, lowest
                    plotting(task,fped,cped) = out.threshold;
                catch
                    plotting(task,fped,cped) = -1;
                end
            end
        end
    end
    hold off
    figure
    subplot(211)
    imagesc(squeeze(plotting(1,:,:)));
    colorbar
    xlabel('Flow Pedestals');
    ylabel('Contrast Pedestals');
    ax = gca;
    ax.XTick = 1:4;
    ax.XTickLabel = stimulus.pedestals.flow;
    ax.YTick = 1:4;
    ax.YTickLabel = stimulus.pedestals.contrast;
    title('Task: Flow')
    subplot(212)
    imagesc(squeeze(plotting(2,:,:)));
    colorbar
    xlabel('Flow Pedestals');
    ylabel('Contrast Pedestals');
    ax = gca;
    ax.XTick = 1:4;
    ax.XTickLabel = stimulus.pedestals.flow;
    ax.YTick = 1:4;
    ax.YTickLabel = stimulus.pedestals.contrast;
    title('Task: Contrast')
%     hold on
%     title(sprintf('%s, R->G->B High',taskOpts{task}));
%     plot(stimulus.pedestals.(taskOpts{task})(1:4),plotting(1,:),'-r');
%     plot(stimulus.pedestals.(taskOpts{task})(1:4),plotting(2,:),'--r');
%     legend(taskOpts);
%     axis([stimulus.pedestals.(taskOpts{task})(1) stimulus.pedestals.(taskOpts{task})(4) 0 .3]);
%     hold off

catch
    disp('(flowAwe) Figures were not generated successfully.');
end

%% checkStaircaseStop
function checkStaircaseStop(stimulus)
% Check both staircases
for task = 1:2
    for ped = 1:4
        for ped2 = 1:4
            s = stimulus.staircase{task,ped,ped2};
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
                stimulus.staircase{task,ped,ped2} = doStaircase('init',s,'initialThreshold',vals{threshPos},'initialStepsize',vals{stepPos});
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create dots for optic flow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = initDotsOpticflow(dots,myscreen)

% focal length to projection plane
% projection plane is defined to be 
% 1 unit wide and high, so with 
% this focal length, we are looking at
% a view of the world with a 90 deg fov
dots.f = .5;

% translation and rotation matrices
dots.T = [0 0 dots.speed/myscreen.framesPerSecond];
dots.R = [0 0 0];

% maximum depth of points
dots.maxZ = 20;dots.minZ = dots.f;
dots.maxX = 20;
dots.maxY = 20;

% make a brick of points
dots.n = round(myscreen.imageWidth*myscreen.imageHeight*dots.density);
% make sure it's an even number
dots.n = dots.n + mod(dots.n,2);

% set half to white and half to black
dots.con = repmat([1 2],1,dots.n/2);

% initial position of dots
dots.X = 25 + dots.mult*abs(2*dots.maxX*rand(1,dots.n)-dots.maxX);
dots.Y = 2*dots.maxY*rand(1,dots.n)-dots.maxY;
dots.Z = (dots.maxZ-dots.minZ)*rand(1,dots.n)+dots.minZ;

% get projection on to plane
dots.xproj = dots.f*dots.X./dots.Z;
dots.yproj = dots.f*dots.Y./dots.Z;

% put into screen coordinates
dots.x = dots.xproj*myscreen.imageWidth;
dots.y = dots.yproj*myscreen.imageHeight;

% set incoherent dots to 0
dots.coherency = 1;
dots.incoherent = rand(1,dots.n) > dots.coherency;
dots.incoherentn = sum(dots.incoherent);
dots.coherent = ~dots.incoherent;

dots.randT = zeros(3,dots.incoherentn);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% step dots for opticflow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dots = updateDotsOpticFlow(dots,coherence,myscreen)

% get the coherent and incoherent dots
% if (dots.coherency ~= coherence)
  dots.incoherent = rand(1,dots.n) > coherence;
  dots.incoherentn = sum(dots.incoherent);
  dots.coherent = ~dots.incoherent;
  dots.coherency = coherence;
  % generate a random transformation matrix for each incoherent point
  dots.randT = rand(3,dots.incoherentn)-0.5;
  % and normalize the transformation to have the same length
  % (i.e. speed) as the real transformation matrix
  dots.randT = sqrt(sum(dots.T.^2))*dots.randT./([1 1 1]'*sqrt(sum(dots.randT.^2)));
% end

% update relative position of dots in 3-space to observer
dots.X(dots.coherent) = dots.X(dots.coherent)-dots.T(1)*dots.dir;
dots.Y(dots.coherent) = dots.Y(dots.coherent)-dots.T(2)*dots.dir;
dots.Z(dots.coherent) = dots.Z(dots.coherent)-dots.T(3)*dots.dir*dots.speedmult;

% now move the incoherent points according to the random trasnformation
dots.X(dots.incoherent) = dots.X(dots.incoherent)-dots.randT(1,:);
dots.Y(dots.incoherent) = dots.Y(dots.incoherent)-dots.randT(2,:);
dots.Z(dots.incoherent) = dots.Z(dots.incoherent)-dots.randT(3,:);

% get all points that have fallen off the screen
offscreen = dots.Z<dots.minZ;

% and put them at the furthest distance
dots.Z(offscreen) = dots.maxZ;

% get all points that have fallen out of view
offscreen = dots.Z>dots.maxZ;
% and move them to the front plane
dots.Z(offscreen) = dots.minZ;

% put points fallen off the X edge back
if dots.mult > 0 % we are looking at RIGHT dots
    offscreen = dots.X < dots.centerOffset * dots.mult;
    dots.X(offscreen) = dots.X(offscreen)+dots.maxX;
    offscreen = dots.X > dots.maxX + dots.centerOffset*dots.mult;
    dots.X(offscreen) = dots.X(offscreen)-dots.maxX;
else % LEFT dots
    offscreen = dots.X < -dots.maxX + dots.centerOffset * dots.mult;
    dots.X(offscreen) = dots.X(offscreen)+dots.maxX;
    offscreen = dots.X > dots.centerOffset * dots.mult;
    dots.X(offscreen) = dots.X(offscreen)-dots.maxX;
end

% put points fallen off the Y edge back
offscreen = dots.Y < -dots.maxY;
dots.Y(offscreen) = dots.Y(offscreen)+2*dots.maxY;
offscreen = dots.Y > dots.maxY;
dots.Y(offscreen) = dots.Y(offscreen)-2*dots.maxY;

% project on to screen
dots.xproj = dots.f*dots.X./dots.Z;
dots.yproj = dots.f*dots.Y./dots.Z;

% stuff to compute median speed
dots.oldx = dots.x;
dots.oldy = dots.y;

% get actual screen coordinates
dots.x = dots.xproj*myscreen.imageWidth;
dots.y = dots.yproj*myscreen.imageHeight;