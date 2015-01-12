% cue_gender
%
%      usage: Called through cue_noisecon
%         by: daniel birman
%       date: 11/10/14
%    purpose: gender discrimination of peripheral stimuli, relies on flags
%    from cue_noisecon to display.
%

function [task, myscreen] = cue_gender(myscreen)
%% Initialize Stimulus

% NOTE: MALE = 8, FEMALE = 9. Note that this depends on the order of the
% folders.
flip = 0; %when flipped, 9 = FEMALE

global stimulus
stimulus.p.lastpick = zeros(1,4);
stimulus.p.posx = [10 -10];
stimulus.p.posy = [0 0];
stimulus.p.numImages(1) = length(stimulus.p.tex{1});
stimulus.p.numImages(2) = length(stimulus.p.tex{2});
if flip
    stimulus.p.responseLetters = [8 9 10];
else
    stimulus.p.responseLetters = [9 8 10];
end
stimulus.p.init_SOA = .25;
stimulus.p.scram.rate = 12.5;
stimulus.p.scram.last = 0;
stimulus.p.scram.left = 0;
stimulus.p.scram.right = 0;

%% Task Params

% Segments are: NOTHING, MASK STREAM, PRESENTATION, PRESENTATION MASK
task{1}.segmin = [inf .010 inf 1]; 
task{1}.segmax = [inf .225 inf 1]; 
% We only get responses after presentation
task{1}.getResponse = [0 0 1 1];
task{1}.randVars.calculated.position = nan;
task{1}.randVars.calculated.intervals = nan; 
task{1}.randVars.calculated.gender = nan(1,2);
task{1}.randVars.calculated.images = nan(1,2);
task{1}.randVars.calculated.respond = nan;
task{1}.randVars.calculated.mainTrialNum = nan;
task{1}.randVars.calculated.SOA = nan;
task{1}.randVars.calculated.sOnset = nan;
task{1}.waitForBacktick = 1;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% init staircase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if stimulus.initStair
    threshold = stimulus.p.init_SOA;
    useLevittRule = 1;
    disp(sprintf('(noisecon) Initializing staircase with threshold: %f stepsize: %f useLevittRule: %i',threshold,useLevittRule));
    stimulus = initStaircase(threshold,stimulus);
end


%% Initialize Task
[task{1}, myscreen] = initTask(task{1},myscreen,@startSegmentCallback,@screenUpdateCallback,@responseCallback,@trialCallback);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Segment Callback %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = startSegmentCallback(task, myscreen)
global stimulus

if task.thistrial.thisseg == 2
    % Everything is random to get a nice scramble, no need to track.
    stimulus.p.g1 = task.thistrial.gender(1);
    stimulus.p.g2 = task.thistrial.gender(2);
    
    stimulus.p.n1 = task.thistrial.images(1);
    stimulus.p.n2 = task.thistrial.images(2);
    
    stimulus.p.scramble = 1;
elseif task.thistrial.thisseg == 3
    if task.thistrial.display == 1
        % We should do this interval
        task.thistrial.respond = 1;
        stimulus.p.SOA_onset{task.trialnum} = mglGetSecs;
        stimulus.p.scramble = 0;
    else
        stimulus.p.scramble = 1;
        stimulus.p.SOA_onset{task.trialnum} = mglGetSecs;
        task = jumpSegment(task);
    end
else
    stimulus.p.scramble = 1;
end

% **************************IN GET RESPONSE CALLBACK**********************************
function [task, myscreen] = responseCallback(task, myscreen)
global stimulus

if any(task.thistrial.whichButton == stimulus.p.responseLetters)
    if task.thistrial.gotResponse == 0 
        whichGender = task.thistrial.gender(task.thistrial.position);
        if task.thistrial.respond == 0
            % Nothing displayed
            if task.thistrial.whichButton == stimulus.p.responseLetters(3)
                % Participant says they saw nothing (correct)                
                correctIncorrect = 'correct';
                if stimulus.dual
                    stimulus.p.dualstaircase{stimulus.blocks.curBlock} = doStaircase('update',stimulus.p.dualstaircase{stimulus.blocks.curBlock},1);
                else
                    stimulus.p.staircase = doStaircase('update',stimulus.p.staircase,1);
                end
                mglFillOval(0, 0, [.5 .5],  stimulus.colors.reservedColor(15));
            else
                correctIncorrect = 'incorrect';
                % Participant saw something (incorrect)
                if stimulus.dual
                    stimulus.p.dualstaircase{stimulus.blocks.curBlock} = doStaircase('update',stimulus.p.dualstaircase{stimulus.blocks.curBlock},0);
                else
                    stimulus.p.staircase = doStaircase('update',stimulus.p.staircase,0);
                end
                mglFillOval(0, 0, [.5 .5],  stimulus.colors.reservedColor(14));
            end
        else
            % Something displayed, check if they got it correct
            if (task.thistrial.whichButton == stimulus.p.responseLetters(whichGender))
                correctIncorrect = 'correct';
                if stimulus.dual
                    stimulus.p.dualstaircase{stimulus.blocks.curBlock} = doStaircase('update',stimulus.p.dualstaircase{stimulus.blocks.curBlock},1);
                else
                    stimulus.p.staircase = doStaircase('update',stimulus.p.staircase,1);
                end
                mglFillOval(0, 0, [.5 .5],  stimulus.colors.reservedColor(15));
            else
                correctIncorrect = 'incorrect';
                if stimulus.dual
                    stimulus.p.dualstaircase{stimulus.blocks.curBlock} = doStaircase('update',stimulus.p.dualstaircase{stimulus.blocks.curBlock},0);
                else
                    stimulus.p.staircase = doStaircase('update',stimulus.p.staircase,0);
                end
                mglFillOval(0, 0, [.5 .5],  stimulus.colors.reservedColor(14));
            end
        end
        disp(sprintf('(gender) Response %s',correctIncorrect));
    else
        disp(sprintf('(gender) Multiple responses... (%i ignored)',task.thistrial.whichButton));
    end
end

% ********************************screen update callback  ******************************

function [task, myscreen] = screenUpdateCallback(task, myscreen)
global stimulus

% at every screen refresh, check the flags:
if stimulus.pFlag == 2
    % Do nothing
    return
elseif stimulus.pFlag == 1  % If the main task has set the start flag to 1
  stimulus.pFlag = 0;   % reset it to 0
  task = jumpSegment(task);       %  and start the subsidiary task by jumping to the next segment
end

if task.thistrial.thisseg == 1
    return
end

if stimulus.p.scram.left == 0 || (mglGetSecs - stimulus.p.scram.last) > (1 / stimulus.p.scram.rate)
    stimulus.p.scram.last = mglGetSecs;
    old = [stimulus.p.scram.left stimulus.p.scram.right];
    while any([stimulus.p.scram.left == old, stimulus.p.scram.right == old])
        stimulus.p.scram.left = randi(9) + 1;
        stimulus.p.scram.right = randi(9) + 1;
    end
end
    
if stimulus.p.scramble == 1 
% Just draw random stuff
    imgDraw(stimulus.p.g1,stimulus.p.n1,1,1,stimulus.p.scram.left,stimulus);
    imgDraw(stimulus.p.g2,stimulus.p.n2,2,1,stimulus.p.scram.right,stimulus);
else
    % Draw the correct image on one side, scrambled on the other
    imgDraw(stimulus.p.g1,stimulus.p.n1,1,task.thistrial.position==1,stimulus.p.scram.left,stimulus);
    imgDraw(stimulus.p.g2,stimulus.p.n2,2,task.thistrial.position==2,stimulus.p.scram.right,stimulus);
end

% Jump segment when the timing calls for it
if task.thistrial.thisseg == 3 && (mglGetSecs - stimulus.p.SOA_onset{task.trialnum}) > task.thistrial.SOA
    task = jumpSegment(task);
end

function imgDraw(gen,imgN,pos,scramble,s,stimulus)

if scramble == 1
    % get a mask image
    image = stimulus.p.tex{gen,s}(imgN);
else
    % actual image
    image = stimulus.p.tex{gen,1}(imgN);
end

% Move to the correct position
mglBltTexture(image,[stimulus.p.posx(pos) stimulus.p.posy(pos) stimulus.p.widthDeg stimulus.p.heightDeg]);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Trial Callback %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = trialCallback(task, myscreen)
global stimulus

task.thistrial.mainTrialNum = stimulus.curTrial;
    
stimulus.p.scram.left = 0;
stimulus.p.scram.right = 0;

% by default, don't look for responses later
task.thistrial.respond = 0;
% get the SOA for this trial
if stimulus.dual
    [task.thistrial.SOA, stimulus.p.dualstaircase{stimulus.blocks.curBlock}] = doStaircase('testValue',stimulus.p.dualstaircase{stimulus.blocks.curBlock});
else
    [task.thistrial.SOA, stimulus.p.staircase] = doStaircase('testValue',stimulus.p.staircase);
end

genBase = [1 2];
task.thistrial.position = randi(2);
intOpts = [0 1 1 1 1];
task.thistrial.display = intOpts(randi(length(intOpts)));
task.thistrial.gender = genBase(randperm(2));
task.thistrial.images(1) = randi(stimulus.p.numImages(task.thistrial.gender(1)));
task.thistrial.images(2) = randi(stimulus.p.numImages(task.thistrial.gender(2)));

%%%%%%%%%%%%%%%%%%%%%%%%
%    initStaircase     %
%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initStaircase(threshold,stimulus)

stimulus.p.staircase = doStaircase('init','upDown','initialThreshold', ...
    threshold,'initialStepsize',threshold/3,'minThreshold=0.01','maxThreshold=.25', ...
    'stepRule','levitt','nTrials=60');
stimulus.p.dualstaircase{1} = stimulus.p.staircase;
stimulus.p.dualstaircase{2} = stimulus.p.staircase;