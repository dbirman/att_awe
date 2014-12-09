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

% NOTE: MALE = 9, FEMALE = 0. Note that this depends on the order of the
% folders.
flip = 0; %when flipped, 9 = FEMALE

global stimulus
stimulus.p.lastpick = zeros(1,4);
stimulus.p.posx = [12 -12];
stimulus.p.posy = [0 0];
stimulus.p.numImages = length(stimulus.p.tex{1});
if flip
    stimulus.p.responseLetters = [9 10];
else
    stimulus.p.responseLetters = [10 9];
end
stimulus.p.init_SOA = .05;   
stimulus.p.scram.rate = 15;
stimulus.p.scram.last = 0;
stimulus.p.scram.left = 0;
stimulus.p.scram.right = 0;

%% Task Params

% Segments are: NOTHING, MASK STREAM, PRESENTATION, PRESENTATION MASK
task{1}.segmin = [inf .010 inf 1]; 
task{1}.segmax = [inf .3 inf 1]; 
% We only get responses after presentation
task{1}.getResponse = [0 0 1 1];
task{1}.randVars.calculated.position = nan;
task{1}.randVars.calculated.intervals = nan; % Note, you can set this to 0/1/2 to have the images only show up during some of the intervals
task{1}.randVars.calculated.gender = nan(1,2);
task{1}.randVars.calculated.images = nan(1,2);
task{1}.randVars.calculated.respond = nan;
task{1}.randVars.calculated.SOA = nan;
task{1}.randVars.calculated.sOnset = nan;
task{1}.waitForBacktick = 1;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% init staircase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if stimulus.initStair
    threshold = stimulus.p.init_SOA;
    stepsize = .005;
    useLevittRule = 1;
    disp(sprintf('(noisecon) Initializing staircase with threshold: %f stepsize: %f useLevittRule: %i',threshold,stepsize,useLevittRule));
    stimulus = initStaircase(threshold,stimulus,stepsize,useLevittRule);
end


%% Initialize Task
[task{1}, myscreen] = initTask(task{1},myscreen,@startSegmentCallback,@screenUpdateCallback,@responseCallback,@trialCallback);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Segment Callback %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = startSegmentCallback(task, myscreen)
global stimulus

curGender = task.thistrial.gender(task.thistrial.position);
if task.thistrial.thisseg == 2
    task.thistrial.SOA = stimulus.staircase{curGender}.threshold;
    % Everything is random to get a nice scramble, no need to track.
    stimulus.p.g1 = randi(2);
    stimulus.p.g2 = randi(2);
    stimulus.p.n1 = randi(stimulus.p.numImages);
    stimulus.p.n2 = randi(stimulus.p.numImages);
    stimulus.p.scramble = 1;
elseif task.thistrial.thisseg == 3
    % Tracking is important
    stimulus.p.g1 = task.thistrial.gender(1); % gender(POSITION)
    stimulus.p.g2 = task.thistrial.gender(2);
    
    stimulus.p.cur_SOA = stimulus.p.staircase{curGender}.threshold;
    
    stimulus.p.n1 = task.thistrial.images(1);
    stimulus.p.n2 = task.thistrial.images(2);
    
    % Presentation (Choose an image)
% % % % % % %     if task.thistrial.intervals == 1 && stimulus.pInt == 1
% % % % % % %         % There is only one interval and we are in the first one
% % % % % % %         task.thistrial.respond = 1;
% % % % % % %         stimulus.p.SOA_onset{task.trialnum}(1) = mglGetSecs;
% % % % % % %         stimulus.p.scramble = 0;
% % % % % % %         % One interval, second one
% % % % % % %     elseif task.thistrial.intervals == 2 && stimulus.pInt == 2
% % % % % % %         task.thistrial.respond = 1;
% % % % % % %         stimulus.p.SOA_onset{task.trialnum}(2) = mglGetSecs;
% % % % % % %         stimulus.p.scramble = 0;
% % % % % % %     elseif task.thistrial.intervals == 3
    if task.thistrial.intervals == 3
        % We should do both intervals
        task.thistrial.respond = 1;
        stimulus.p.SOA_onset{task.trialnum} = mglGetSecs;
        stimulus.p.scramble = 0;
    else
%         % 0 response trials, we don't care which we are in.
%         % We just skip to the mask phase immediately
%         task.thistrial.respond = 0;
%         stimulus.p.SOA_onset{task.trialnum} = mglGetSecs;
%         task = jumpSegment(task);
        warning('Shouldn''t get here');
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
        if (task.thistrial.whichButton == stimulus.p.responseLetters(whichGender))
            correctIncorrect = 'correct';
            stimulus.p.staircase{whichGender} = upDownStaircase(stimulus.p.staircase{whichGender},1);
            mglFillOval(0, 0, [.5 .5],  stimulus.colors.reservedColor(15));
        else
            correctIncorrect = 'incorrect';
            stimulus.p.staircase{whichGender} = upDownStaircase(stimulus.p.staircase{whichGender},0);
            mglFillOval(0, 0, [.5 .5],  stimulus.colors.reservedColor(14));
        end
        disp(sprintf('(Peripheral) Response %s',correctIncorrect));
    else
        disp(sprintf('(Peripheral) Multiple responses... (%i ignored)',task.thistrial.whichButton));
    end
end

% ********************************screen update callback  ******************************

function [task, myscreen] = screenUpdateCallback(task, myscreen)
global stimulus

% % mglTextDraw(num2str(task.thistrial.thisseg),[0,-8]);

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
if task.thistrial.thisseg == 3 && (mglGetSecs - stimulus.p.SOA_onset{task.trialnum}) > stimulus.p.cur_SOA
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
mglBltTexture(image,[stimulus.p.posx(pos) stimulus.p.posy(pos) stimulus.widthDeg stimulus.heightDeg]);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Trial Callback %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = trialCallback(task, myscreen)
global stimulus

% by default, don't look for responses later
task.thistrial.respond = 0;

genBase = [1 2];
task.thistrial.position = randi(2);
task.thistrial.intervals = 3;
task.thistrial.gender = genBase(randperm(2));
task.thistrial.images = randi(stimulus.p.numImages,1,2);

%%%%%%%%%%%%%%%%%%%%%%%%
%    initStaircase     %
%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initStaircase(threshold,stimulus,stepsize,useLevittRule)

for gender = 1:2 % female / male
    stimulus.p.staircase{gender} = upDownStaircase(1,2,threshold,stepsize,useLevittRule);
    stimulus.p.staircase{gender}.minThreshold = 0;
    stimulus.p.staircase{gender}.maxThreshold = .2;
end