% flowAwe
%
%      usage: myscreen=flowAwe()
%         by: daniel birman
%       date: 11/10/14
%    purpose: contrast change detection with cued selective attention.
%
%        use: call cue_fade() to initialize. The first
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

function [myscreen] = cue_fade(varargin)

global stimulus
%% Initialize Variables

% add arguments later
widthDeg = []; heightDeg = [];
stimFileNum = [];
plots = [];
practice = [];
getArgs(varargin,{'widthDeg=5.5', 'heightDeg=5.5', ...
    'stimFileNum=-1', ...
    'dual=0','plots=1','practice=0'});

stimulus.counter = 1; % This keeps track of what "run" we are on.
%% Setup Screen

screen.screenNumber = 2;
myscreen = initScreen(screen);

%% Open Old Stimfile
stimulus.initStair = 1;

if ~isempty(mglGetSID) && isdir(sprintf('~/data/cue_fade/%s',mglGetSID))
    % Directory exists, check for a stimefile
    files = dir(sprintf('~/data/cue_fade/%s/1*mat',mglGetSID));

    if length(files) >= 1
        if stimFileNum == -1
            if length(files) > 1
                warning('Multiple stimfiles found, loading last one. If you wanted a different functionality use stimFileNum=#');
            end
            fname = files(end).name;
        else
            fname = files(stimFileNum).name;
        end
        s = load(sprintf('~/data/cue_fade/%s/%s',mglGetSID,fname));
        stimulus.staircase = s.stimulus.staircase;
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

stimulus.colors.rmed = 127.75;
stimulus.colors.mrmax = 255;
stimulus.colors.mrmin = 0;

stimulus.wedge.deg = 15;

stimulus.dots.xcenter = 0;
stimulus.dots.ycenter = 0;
stimulus.dots.dotsize = 2;
stimulus.dots.density = 8;
stimulus.dots.coherence = 1;
stimulus.dots.speed = 4;
stimulus.dots.T = [0 0 stimulus.dots.speed/myscreen.framesPerSecond];
stimulus.dots.dir = 0;
stimulus.dotsR = stimulus.dots;
stimulus.dotsR.mult = 1;
stimulus.dotsL = stimulus.dots;
stimulus.dotsL.mult = -1;
stimulus = rmfield(stimulus,'dots');

stimulus.mask = 1;

stimulus.pedestals.pedOpts = {'flow','contrast'};
stimulus.pedestals.flow = [.25 .5 .75];
stimulus.pedestals.initThresh.flow = .2;
stimulus.pedestals.contrast = [.25 .5 .75]; % for now...
stimulus.pedestals.initThresh.contrast = .2; 

stimulus.dotsR = initDotsOpticflow(stimulus.dotsR,myscreen);
stimulus.dotsL = initDotsOpticflow(stimulus.dotsL,myscreen);

%% MGL Parameters
mglTextSet('Helvetica',32,1,0,0,0,0,0,0,0);
    
%% Additional Dots Params

% set color
stimulus.dotsR.color = ones(stimulus.dotsR.n,1);
stimulus.dotsL.color = ones(stimulus.dotsL.n,1);

% create stencil
if stimulus.mask
  mglClearScreen;
  mglStencilCreateBegin(1,1);
  % and draw that oval
  % now draw the wedge
  mglGluPartialDisk(0,0,0,20,-stimulus.wedge.deg,2*stimulus.wedge.deg,[1 1 1],16);
  mglGluPartialDisk(0,0,0,20,180-stimulus.wedge.deg,2*stimulus.wedge.deg,[1 1 1],16);
  mglGluPartialDisk(0,0,0,2,0,360,[1 1 1],60);
  mglStencilCreateEnd;
  mglClearScreen;
end

%% Setup Task

% This is the contrast change detection task
task{1}{1}.waitForBacktick = 1;
stimulus.seg.ITI = 1; % the ITI is either 20s (first time) or 1s
stimulus.seg.stim = 2;
stimulus.seg.ISI = 3;
stimulus.seg.resp = 4;
task{1}{1}.segmin = [1 1 .1 1.2];
task{1}{1}.segmax = [1 1 .5 1.2];
task{1}{1}.synchToVol = [0 0 0 0];
task{1}{1}.getResponse = [0 0 0 1];
task{1}{1}.parameter.side = [1 2]; % 1 = left, 2 = right, the side will be the one with con/flow + delta (From staircase)
task{1}{1}.parameter.conPedestal = [1 2 3]; % target contrast
task{1}{1}.parameter.floPedestal = [1 2 3]; % target flow coherence
task{1}{1}.random = 1;
task{1}{1}.numBlocks = 6;

%% Run variables
task{1}{1}.randVars.calculated.task = nan; % Current task (calc per run)
task{1}{1}.randVars.calculated.deltaPed = nan; % Current task (calc per run)
task{1}{1}.randVars.calculated.coherence = nan;
task{1}{1}.randVars.calculated.contrast = nan;

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
stimulus.runs.curTask = stimulus.runs.taskList(stimulus.counter);

%% Full Setup
% Initialize task (note phase == 1)
for phaseNum = 1:length(task{1})
    [task{1}{phaseNum}, myscreen] = initTask(task{1}{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@getResponseCallback,@startTrialCallback,[],@startBlockCallback);
end

mglClearScreen(stimulus.colors.rmed/stimulus.colors.mrmax);

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
mglClearScreen(stimulus.colors.rmed/stimulus.colors.mrmax);
mglTextDraw(stimulus.runs.taskOptsText{stimulus.runs.curTask},[0 0]);
mglFlush

myscreen.flushMode = 1;

% let the user know
disp(sprintf('(flowAwe) Starting run number: %i',stimulus.counter));

%% Main Task Loop

phaseNum = 1;
% Again, only one phase.
while (phaseNum <= length(task{1})) && ~myscreen.userHitEsc
    % update the task
    [task{1}, myscreen, phaseNum] = updateTask(task{1},myscreen,phaseNum);
    % flip screen
    myscreen = tickScreen(myscreen,task);
end

% task ended
mglClearScreen(stimulus.colors.rmed/stimulus.colors.mrmax);
mglTextDraw('Run complete... please wait.',[0 0]);
mglFlush
myscreen.flushMode = 1;
disp('(flowAwe) Run ending...');

if plots
    disp('(flowAwe) Displaying plots');
    dispStaircase(stimulus);
end

% delete texture
% if isfield(stimulus,'tex')
%     for x = 1:size(stimulus.tex,1)
%         for y = 1:size(stimulus.tex,2)
%             mglDeleteTexture(stimulus.tex(x,y));
%         end
%     end
% end
% mglDeleteTexture(stimulus.mask);
% stimulus = rmfield(stimulus,'tex');

%% Save file (first file will always have full textures)
if stimulus.counter > 1
    % Temporarily remove raw files to limit stimfile size, then replace these
    % so they can be used on the next run.
    stimbackup = stimulus;
    stimulus = rmfield(stimulus,'raw');
    stimulus = rmfield(stimulus,'images');
end

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);

if stimulus.counter > 1
    stimulus = stimbackup;
    clear stimbackup
end

%%%%%%%%%%%%%%%%%%%%%%%%% EXPERIMENT OVER: HELPER FUNCTIONS FOLLOW %%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Runs at the start of each Trial %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = startTrialCallback(task,myscreen)

myscreen.flushMode = 0;
% We need to do a number of things here. First we need to set up the
% missing thistrial parameters. Second, we will pre-calculate the textures
% that will be displayed on this particular trial. 

global stimulus

task.thistrial.coherence = stimulus.pedestals.flow(task.thistrial.floPedestal);
task.thistrial.contrast = stimulus.pedestals.contrast(task.thistrial.conPedestal);
[task.thistrial.deltaPed, stimulus] = getDeltaPed(stimulus,stimulus.runs.curTask,curPedValue(task,stimulus));

if stimulus.runs.curTask==1
    % coherence
    stimulus.live.cohDelta = task.thistrial.deltaPed;
    stimulus.live.conDelta = 0;
    disp(sprintf('(flowAwe) Trial ## starting. Coherence: %.02f + %.02f Contrast %.02f',task.thistrial.coherence,stimulus.live.cohDelta,task.thistrial.contrast));
else
    % contrast
    stimulus.live.cohDelta = 0;
    stimulus.live.conDelta = task.thistrial.deltaPed;
    disp(sprintf('(flowAwe) Trial ## starting. Coherence: %.02f Contrast %.02f + $.02f',task.thistrial.coherence,task.thistrial.contrast,stimulus.live.conDelta));
end

function ped = curPedValue(task,stimulus)
if stimulus.runs.curTask==1
    ped = task.thistrial.floPedestal;
else
    ped = task.thistrial.conPedestal;
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Runs at the start of each Segment %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = startSegmentCallback(task, myscreen)

global stimulus

switch task.thistrial.thisseg
    case stimulus.seg.ITI
        stimulus.live.dots = 0;
        stimulus.live.fixColor = 0;
    case stimulus.seg.stim
        stimulus.live.dots = 1;
        stimulus.live.fixColor = 0;
    case stimulus.seg.ISI
        stimulus.live.dots = 0;
        stimulus.live.fixColor = 0;
    case stimulus.seg.resp
        stimulus.live.dots = 0;
        stimulus.live.fixColor = 1;
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

mglClearScreen(stimulus.colors.rmed/stimulus.colors.mrmax);

upFix(stimulus);
if stimulus.live.dots==1, stimulus = upDots(task,stimulus,myscreen); end;

%%
function upFix(stimulus)

mglFixationCross(1,1,stimulus.live.fixColor);

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

% draw the dots
if stimulus.mask,mglStencilSelect(1);end
% dotsR
% update +contrast
mglPoints2(stimulus.dotsR.x(stimulus.dotsR.con==1),stimulus.dotsR.y(stimulus.dotsR.con==1),...
    stimulus.dotsR.dotsize,[.5 .5 .5] - task.thistrial.contrast/2 - rConDelta);
% update - contrast
mglPoints2(stimulus.dotsR.x(stimulus.dotsR.con==2),stimulus.dotsR.y(stimulus.dotsR.con==2),...
    stimulus.dotsR.dotsize,[.5 .5 .5] + task.thistrial.contrast/2 + rConDelta);
% dotsL
% update +contrast
mglPoints2(stimulus.dotsL.x(stimulus.dotsL.con==1),stimulus.dotsL.y(stimulus.dotsL.con==1),...
    stimulus.dotsL.dotsize,[.5 .5 .5] - task.thistrial.contrast/2 - lConDelta);
% update - contrast
mglPoints2(stimulus.dotsL.x(stimulus.dotsL.con==2),stimulus.dotsL.y(stimulus.dotsL.con==2),...
    stimulus.dotsL.dotsize,[.5 .5 .5] + task.thistrial.contrast/2 + lConDelta);
if stimulus.mask,mglStencilSelect(0);end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Called When a Response Occurs %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function [task, myscreen] = getResponseCallback(task, myscreen)

global stimulus

responseText = {'Incorrect','Correct'};
fixColors = {[1 0 0],[0 1 0]};

if any(task.thistrial.whichButton == stimulus.responseKeys)
    if task.thistrial.gotResponse == 0
        task.thistrial.correct = task.thistrial.whichButton == stimulus.responseKeys(task.thistrial.side);
        % Store whether this was correct
        stimulus.live.fixColor = fixColors{task.thistrial.correct+1};
        disp(sprintf('(flowAwe) Response %s',responseText{task.thistrial.correct+1}));
%         if ~task.thistrial.catch
            stimulus.staircase{stimulus.runs.curTask,curPedValue(task,stimulus)} = ...
                doStaircase('update',stimulus.staircase{stimulus.runs.curTask,curPedValue(task,stimulus)},task.thistrial.correct);
%         else
    else
        disp(sprintf('(flowAwe) Subject responded multiple times: %i',task.thistrial.gotResponse+1));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Block Call Back %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function [task, myscreen] = startBlockCallback(task, myscreen)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                              HELPER FUNCTIONS                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% getDeltaPed

function [deltaPed, stimulus] = getDeltaPed(stimulus,taskNum,pedNum)

[deltaPed, stimulus.staircase{taskNum,pedNum}] = doStaircase('testValue',stimulus.staircase{taskNum,pedNum});


%%%%%%%%%%%%%%%%%%%%%%%%
%    initStaircase     %
%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initStaircase(stimulus)

stimulus.staircase = cell(2,length(stimulus.pedestals.contrast));
for task = 1:2
    stimulus.staircase{task,1} = doStaircase('init','upDown',...
        'initialThreshold',stimulus.pedestals.initThresh.(stimulus.pedestals.pedOpts{task}),...
        'initialStepsize',stimulus.pedestals.initThresh.(stimulus.pedestals.pedOpts{task})/3,...
        'minThreshold=0.001','maxThreshold=0.5','stepRule','pest', ...
        'nTrials=50','maxStepsize=.2','minStepsize=.001');
end

for i = 2:length(stimulus.pedestals.contrast)
    stimulus.staircase{1,i} = stimulus.staircase{1,1};
    stimulus.staircase{2,i} = stimulus.staircase{2,1};
end

%%%%%%%%%%%%%%%%%%%%%%%
%    dispStaircase    %
%%%%%%%%%%%%%%%%%%%%%%%
function dispStaircase(stimulus)

try
    taskOpts = {'flow','contrast'};
    
    figure % this is the 'staircase' figure
    title('%s, Staircase plot (R->G->B high)');
    hold on
    drawing = {'-r' '-g' '-b'
                '--r' '--g' '--b'};
    for task = 1:2
        for ped = 1:3
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
    title(sprintf('%s, R->G->B High'));
    plot(stimulus.pedestals.(taskOpts{task})(1:3),plotting(1,:),'-r');
    plot(stimulus.pedestals.(taskOpts{task})(1:3),plotting(2,:),'--r');
    axis([stimulus.pedestals.(taskOpts{task})(1) stimulus.pedestals.(taskOpts{task})(3) 0 1]);
    hold off

catch
    disp('(flowAwe) Figures were not generated successfully.');
end

%% checkStaircaseStop
function checkStaircaseStop(stimulus)
% Check both staircases
for cues = 1:2
    for ped = 1:3
        s = stimulus.staircase{cues,ped,stimulus.runs.curTask};
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
            disp('THIS CODE IS NOT CERTAIN TO WORK! CHECK THE OUTPUT!');
            stimulus.staircase{cues,ped,stimulus.runs.curTask} = doStaircase('init',s,'initialThreshold',vals{threshPos},'initialStepsize',vals{stepPos});
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
dots.maxZ = 15;dots.minZ = dots.f;
dots.maxX = 15;
dots.maxY = 15;

% make a brick of points
dots.n = round(myscreen.imageWidth*myscreen.imageHeight*dots.density);

% choose contrast for each point
dots.con = randi(2,1,dots.n);

% initial position of dots
dots.X = dots.mult*abs(2*dots.maxX*rand(1,dots.n)-dots.maxX);
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
dots.X(dots.coherent) = dots.X(dots.coherent)-dots.T(1);
dots.Y(dots.coherent) = dots.Y(dots.coherent)-dots.T(2);
dots.Z(dots.coherent) = dots.Z(dots.coherent)-dots.T(3);

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
offscreen = dots.X < -dots.maxX;
dots.X(offscreen) = dots.X(offscreen)+2*dots.maxX;
offscreen = dots.X > dots.maxX;
dots.X(offscreen) = dots.X(offscreen)-2*dots.maxX;

% put points fallen off the Y edge back
offscreen = dots.Y < -dots.maxY;
dots.Y(offscreen) = dots.Y(offscreen)+2*dots.maxY;
offscreen = dots.Y > dots.maxY;
dots.Y(offscreen) = dots.Y(offscreen)-2*dots.maxY;

% project on to screen
dots.xproj = dots.f*dots.X./dots.Z;
dots.yproj = dots.f*dots.Y./dots.Z;
% dots.ytheta = atan(dots.Y./dots.Z);
% dots.xtheta = atan(dots.X./dots.Z);

% dots.yproj = tan(dots.ytheta)*dots.maxZ;
% dots.xproj = tan(dots.xtheta)*dots.maxZ;

% stuff to compute median speed
dots.oldx = dots.x;
dots.oldy = dots.y;

% get actual screen coordinates
dots.x = dots.xproj*myscreen.imageWidth;
dots.y = dots.yproj*myscreen.imageHeight;