%% Contrast / coherence experiment example
% In this experiment two patches of dots were shown and observers compared
% them to each other. On some trials they compared the _contrast_
% (luminance relative to the background) and on other trials they compared
% the _coherence_ (% of dots moving together).
%
% We measured:
%   (1) BOLD signal response to different levels of contrast/coherence
%   (see: Birman, D., & Gardner, J. L. (2018). A quantitative framework for
%   motion visibility in human cortex. Journal of neurophysiology.)
%   (2) Perceptual sensitivity to differences in contrast/coherence (under
%   review)
%
% We want to link together these two measurements to understand how the
% brain "reads out" signals about contrast and coherence. The script will
% walk you through this.
%
% We've marked locations where you can manipulate the code with *********

%% Load data files

% load the motion visibility framework
% this data, and the accompanying functions are available on the OSF: https://osf.io/s7j9p/
load('data.mat'); 

% load the behavioral data
all = [];
for bi = 1:21
    load(sprintf('data/%i_data.mat',bi));
    all = [all ; adata];
end
adata = all;

%% Part 1: BOLD signal measurements
% What measurements did we make? Let's plot the BOLD signal measurements
% which were made for contrast.

% Pull the BOLD response functions for contrast from one subject

% ***** PARAMETERS ****** %
subject = 1;
roi = 'V1';
% *********************** %

% setup some extra info:
cmap = brewermap(11,'PuOr');
t = 0.25:0.5:40.5;
roiIdx = find(cellfun(@(x) strcmp(x,roi),data{subject}.areas),1);
bold = squeeze(data{subject}.adata{1}.cc.resp(roiIdx,:,:));
con = data{subject}.adata{1}.cc.con;

% Collapse the different coherence values to just the unique contrast
% values
ucon = unique(con);
boldAvg = [];
leg = {};
% Setup the figure
figure(1); clf; hold on
for ui = 1:length(ucon)
    % Compute the average over different coherences
    boldAvg(ui,:) = mean(bold(con==ucon(ui),:));
    % Setup the legend
    leg{ui} = sprintf('25 (baseline) +%2.0f%% contrast',ucon(ui)*100);
    % Plot to the figure (w/ change of color according to contrast)
    p(ui) = plot(t,boldAvg(ui,:),'Color',cmap(5-ui,:));
end

% More figure stuff
title(sprintf('Contrast measurements in %s for subject %i',roi,subject));
xlabel('Time (s)');
ylabel('BOLD signal (%)');
legend(p,leg);

%% Part 2: Perceptual sensitivity
% What measurements did we make? Let's plot the % of times the observer
% said that the "right side" was higher contrast than the left, as a
% function of the actual difference in contrast between the two sides. 

% (Only one subject)

% ***** PARAMETERS ****** %
task = 1; % 2= contrast, 1 = coherence
% *********************** %

%     1       2         3        4      5      6     7      8       9
%   task - basecon - basecoh - conL - conR - cohL - cohR - resp - catch -
%      10      11        12
%   pedcon - pedcoh - correct

% task = 1: coherence, 2: contrast
% resp = 0: left, 1: right

% remove all trials that aren't related to contrast
cdata = adata(adata(:,1)==task,:);
tasks = {'coherence','contrast'};
cmaps = {cmap(8,:),cmap(3,:)};
disp(sprintf('Using ALL subjects we have %i trials for task: %s',size(cdata,1),stask));

stask = tasks{task};
cm = cmaps{task};

% Compute the difference in contrast between the two sides
if task==2
    dcon = cdata(:,5)-cdata(:,4);
else
    dcon = cdata(:,7)-cdata(:,6);
end
% Take quantiles to find bins
qs = quantile(dcon,0:.1:1);
qc = zeros(1,length(qs)-1);
ps = zeros(size(qc));
cis = zeros(size(qc));
% For each quantile, compute the percent of rightward responses
for qi = 1:(length(qs)-1)
    % Tracking stuff
    qmin = qs(qi);
    qmax = qs(qi+1);
    qc(qi) = mean([qmin qmax]);
    % Bin the responses within this bin
    idxs = logical((dcon>qmin).*(dcon<=qmax));
    binResp = cdata(idxs,8);
    if ~isempty(binResp)
        ps(qi) = nanmean(binResp);
        % compute a 95% CI
        ci = bootci(1000,@nanmean,binResp);
        cis(qi) = ci(2)-ps(qi);
    else
        ps(qi) = 0.5;
        cis(qi) = 0;
    end
end

% Plot
figure(2+(task-1)); clf; hold on
plot(qc,ps,'o','MarkerFaceColor',cm,'MarkerEdgeColor','w');
errbar(qc,ps,cis,'-','Color',cm);
axis([-0.4 0.4 0 1]);
xlabel('Difference in contrast (Right - Left, %)');
ylabel('Choices to right (%)');
title(sprintf('Proportion of choices to RIGHT when cued to attend to %s',stask));

%% Part 3: A simple linking model for contrast
% One simple way that the brain might read out information about contrast
% is by comparing the activity of the population of neurons responding to
% each patch of dots (Boynton et al. 1999)

% For each trial, compute the 

%% Part 4: Linking contrast and coherence






%% Some notes:
% This script doesn't include a lot of things that we do in the papers,
% like cross-validation and model comparison. But, I hope that it helps you
% get a better sense of the model and how it was implemented, as well as
% some of the raw data!



%% Helper scripts:

nSIDs = [305 329 43 25 300 346 343 344 338 340 348];
bSIDs = [345 330 337 335 349 354 353 356 334 352]; % behavior only participants
aSIDs = [nSIDs bSIDs];

% copy all the files
for ai = 1:length(aSIDs)
    copyfile(fullfile(datafolder,sprintf('s%03.0f_adata.mat',aSIDs(ai))),fullfile('~/proj/ej_class/data/',sprintf('%i_data.mat',ai)));
end