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

% Now plot the naka-rushton function for this area
figure(2); clf; hold on
x = 0:.01:1;
pred = cohcon_predict(roi,x,0,1);
y = pred.contrastResponse;
plot(x,y,'Color',cmap(3,:));
axis([0 1 0 1.5]);
set(gca,'XTick',[0 0.25 0.5 0.75 1],'XTickLabel',[0 25 50 75 100]);
set(gca,'YTick',[0 0.5 1 1.5]','YTickLabel',{'0%','0.5%','1%','1.5%'});
xlabel('Change in contrast (%)');
ylabel('BOLD Signal (%)');
title('Naka-Rushton function for all subjects');











%% Part 2: Perceptual sensitivity (simulation)
% What measurements did we make? Let's plot the % of times the observer
% said that the "right side" was higher contrast than the left, as a
% function of the actual difference in contrast between the two sides. 

% This part

% (Only one subject)
%     1       2         3        4      5      6     7      8       9
%   task - basecon - basecoh - conL - conR - cohL - cohR - resp - catch -
%      10      11        12
%   pedcon - pedcoh - correct
clear qc ps cis

% setup some stuff
tasks = {'coherence','contrast'};
cmaps = {cmap(8,:),cmap(3,:)};


for task = 1:2
    % remove all trials that aren't related to contrast
    cdata = adata(adata(:,1)==task,:);
    disp(sprintf('Using ALL subjects we have %i trials for task: %s',size(cdata,1),stask));

    % Compute the difference in contrast between the two sides
    dcon = cdata(:,5)-cdata(:,4);
    dcoh = cdata(:,7)-cdata(:,6);

    % Take quantiles to find bins
    qs_con = quantile(dcon,0:.1:1);
    qs_coh = quantile(dcoh,0:.1:1);

    % For each quantile, compute the percent of rightward responses
    for qi = 1:(length(qs)-1)
        % Tracking stuff
        qmin = qs_con(qi);
        qmax = qs_con(qi+1);
        qc{task,2}(qi) = mean([qmin qmax]);
        % Bin the contrast responses within this bin
        conidxs = logical((dcon>qmin).*(dcon<=qmax));

        % Tracking stuff
        qmin = qs_coh(qi);
        qmax = qs_coh(qi+1);
        qc{task,1}(qi) = mean([qmin qmax]);
        % Bin the coherence resopnses
        cohidxs = logical((dcoh>qmin).*(dcoh<=qmax));

        conBinResp = cdata(conidxs,8);
        cohBinResp = cdata(cohidxs,8);

        if ~isempty(binResp)
            % contrast
            ps{task,2}(qi) = nanmean(conBinResp);
            % compute a 95% CI
            ci = bootci(1000,@nanmean,conBinResp);
            cis{task,2}(qi) = ci(2)-ps{task,2}(qi);

            % coherence
            ps{task,1}(qi) = nanmean(cohBinResp);
            % compute a 95% CI
            ci = bootci(1000,@nanmean,cohBinResp);
            cis{task,1}(qi) = ci(2)-ps{task,1}(qi);
        else
            ps{task,2}(qi) = 0.5;
            ps{task,1}(qi) = 0.5;
            cis{task,2}(qi) = 0;
            cis{task,1}(qi) = 0;
        end
    end
end

%% Part 2: Perceptual sensitivity (figures)

% ***** PARAMETERS ****** %
% task 2 = contrast, 1 = coherence
task = 1; 
% *********************** %

stask = tasks{task};
cm = cmaps{task};

figure(3+(task-1)); clf; hold on
plot(qc{task,task},ps{task,task},'o','MarkerFaceColor',cm,'MarkerEdgeColor','w');
errbar(qc{task,task},ps{task,task},cis{task,task},'-','Color',cm);
axis([-0.4 0.4 0 1]);
xlabel(sprintf('Difference in %s (Right - Left, %%)',stask));
ylabel('Choices to right (%)');
title(sprintf('Proportion of choices to RIGHT when responding about %s',stask));




%% Part 3: A simple linking model for contrast
% One simple way that the brain might read out information about contrast
% is by comparing the activity of the population of neurons responding to
% each patch of dots (Boynton et al. 1999)

% Compute the likelihood that a sample from one "neural response" exceeds
% the other, given that the two stimuli are separated by a difference in
% contrast of X. I.e.: phi(0,X,1), where phi is the cumulative normal
% distribution. 

% ***** PARAMETERS ****** %
noise = 0.05;
lapse = 0.1;
% *********************** %

% compute the probability of choosing one contrast over another, for the
% four base contrasts that were shown in the experiment (we'll just average
% these together). Note that the contrasts are always positive relative to
% the base, and that the differences are in the range [0, 0.25];
ubase = [0.325, 0.4, 0.55, 0.85];
x = 0:.001:1.1;
inc = 0:.001:.25;
pred = cohcon_predict('V1',x,0,1);
y = pred.contrastResponse;

clear pred prob
for bi = 1:length(ubase)
    base = ubase(bi);
    
    % compute the response for this base contrast
    resp = base+inc;
    v1_resp = zeros(size(resp));
    for ri = 1:length(resp)
        v1_resp(ri) = y(find(x>=resp(ri),1));
    end
    v1_resp = v1_resp - v1_resp(1);
    prob(bi,:) = normcdf(v1_resp,0,noise);
    prob = (prob-0.5) * (1-lapse) + 0.5;
end

% collapse over the base contrasts (as done in the behavioral data)
prob = mean(prob);


% Plot this onto the contrast figure (fig 4) from before -- you may need to
% re-run the block above to plot that figure.
task = 2
stask = tasks{task};
cm = cmaps{task};

figure(3+(task-1)); clf; hold on
plot(qc{task,task},ps{task,task},'o','MarkerFaceColor',cm,'MarkerEdgeColor','w');
errbar(qc{task,task},ps{task,task},cis{task,task},'-','Color',cm);
axis([-0.4 0.4 0 1]);
xlabel(sprintf('Difference in %s (Right - Left, %%)',stask));
ylabel('Choices to right (%)');
title(sprintf('Proportion of choices to RIGHT when responding about %s',stask));

plot([-fliplr(inc) inc],[1-fliplr(prob) prob],'Color',cm);
legend({'Perceptual sensitivity data','V1 Model'},'Location','nw');

%% Part 4: Linking contrast and coherence
% Flexible readout is one model of how contrast and coherence can be read
% out from sensory representations. The idea is to _weight_ the response in
% different areas according to which task is being performed, then to
% compute the same model as above. (Note that in part 3 the model used only
% the V1 response, and fit this to contrast. Here we introduce a new area,
% MT, to fit the responses to coherence). 

% Before we get to the weighting, we need to know how much the observers
% were influenced by the _irrelevant_ feature during the task. We computed
% this above but didn't plot it. Let's take a look:
figure(5); clf
for task = 1:2
    subplot(1,2,task); hold on
    
    stask = tasks{task};
    cm = cmaps{task};

    plot(qc{task,task},ps{task,task},'o','MarkerFaceColor',cm,'MarkerEdgeColor','w');
    errbar(qc{task,task},ps{task,task},cis{task,task},'-','Color',cm);
    
    flip = [2 1];
    plot(qc{task,flip(task)},ps{task,flip(task)},'o','MarkerFaceColor',cmaps{flip(task)},'MarkerEdgeColor','w');
    errbar(qc{task,flip(task)},ps{task,flip(task)},cis{task,flip(task)},'-','Color',cmaps{flip(task)});
    
    axis([-0.4 0.4 0 1]);
    xlabel('Difference in stimulus strength (Right - Left, %)');
    ylabel('Choices to right (%)');
    title(sprintf('Responding about %s',stask));
end

% Now we can see that actually the model needs to account for (at least)
% two things: (1) the perceptual sensitivity to the feature the observers
% are responding about, but also (2) the insensitivity to the irrelevant
% feature.

% Let's try to fit that model. You now get four additional parameters to
% play with, which control the "weight" placed on V1 and MT in each task.
% Note: the "noise" parameter is relative to the weights (because they both
% scale the response functions in a linear manner) so we don't need to
% explicitly define it any more. 

% ***** PARAMETERS ****** %

% weight in the CONTRAST task:
v1_c = 15;
mt_c = -3;

% weight in the COHERENCE task:
v1_m = -1.5;
mt_m = 20;

lapse = 0.3;
% *********************** %


% compute over each base simultaneously
cbases = [0.325, 0.4, 0.55, 0.85];
mbases = [0.15 0.3 0.45 0.6];

clear v1 mt
x = 0:.001:1.25;
inc = 0:.001:.4;
pred = cohcon_predict({'V1','MT'},x,0,1);
v1(:,2) = pred.contrastResponse(:,1);
mt(:,2) = pred.contrastResponse(:,2);
pred = cohcon_predict({'V1','MT'},0,x,1);
v1(:,1) = pred.coherenceResponse(:,1);
mt(:,1) = pred.coherenceResponse(:,2);

clear pred v1resps mtresps
for bi = 1:length(cbases)
    cbase = cbases(bi);
    mbase = mbases(bi);
    
    % compute responses for contrast/coherence
    cinc = cbase+inc;
    minc = mbase+inc;
    
    clear v1resp mtresp
    v1resp = zeros(length(cinc),2);
    mtresp = v1resp;
    for ri = 1:length(binc)
        % compute contrast responses
        v1resp(ri,2) = v1(find(x>=cinc(ri),1),2);
        mtresp(ri,2) = mt(find(x>=cinc(ri),1),2);
        
        v1resp(ri,1) = v1(find(x>=minc(ri),1),1);
        mtresp(ri,1) = mt(find(x>=minc(ri),1),1);
    end
    v1resp = v1resp-repmat(v1resp(1,:),size(v1resp,1),1);
    mtresp = mtresp-repmat(mtresp(1,:),size(mtresp,1),1);
    v1resps(bi,:,:) = v1resp;
    mtresps(bi,:,:) = mtresp;
end

% for clarity
clear v1resp mtresp
% average over bases
v1resp = squeeze(mean(v1resps));
mtresp = squeeze(mean(mtresps));

% note: we have to do something tricky here. Because we're generating the
% perceptual sensitivity as if the stimuli were uncorrelated we don't want
% to sum the responses together, although in reality the way the model
% works is that each area has only one response (the sum of the contrast
% and coherence response). We keep them separate here so we can calculate
% sensitivity separately according to condition.

% separate the responses according to task

% when reporting contrast, it's "as if" you are only shown differences in 
% contrast:
conresp = [v1resp(:,2) mtresp(:,2)];
% when reporting coherence, it's "as if" you are only shown differences in
% coherence:
cohresp = [v1resp(:,1) mtresp(:,1)];

stims = {cohresp conresp};
reports = {[v1_m;mt_m] [v1_c;mt_c]};

% first dimension in prob is the stimulus 
% second dimension is the report
clear prob
for si = 1:length(stims)
    for ri = 1:length(reports)
        tresp = stims{si}*reports{ri};
        prob(si,ri,:) = normcdf(tresp,0,1);
    end
end

% reduce by lapse rate
prob = (prob-0.5)*(1-lapse/2)+0.5;

% plot by report type
subplot(121);
tp = squeeze(prob(1,1,:))';
plot([-fliplr(inc) inc],[1-fliplr(tp) tp],'Color',cmaps{1});
tp = squeeze(prob(2,1,:))';
plot([-fliplr(inc) inc],[1-fliplr(tp) tp],'Color',cmaps{2});

subplot(122);
tp = squeeze(prob(1,2,:))';
plot([-fliplr(inc) inc],[1-fliplr(tp) tp],'Color',cmaps{1});
tp = squeeze(prob(2,2,:))';
plot([-fliplr(inc) inc],[1-fliplr(tp) tp],'Color',cmaps{2});

%% Some notes:
% This script doesn't include a lot of things that we do in the papers!
% 
% For example:
%   To fit the model parameters via maximum-likelihood, how would you 
%   compute the likelihood of a trial?
%   
%   We use cross-validation to ensure that a model is not overfitting. What
%   would that entail here?
%
%   You might want to add various kinds of bias, for example participants
%   might just be more likely to say "left" or "right" on a trial. Do you
%   see any of this in the average data? How could you add it to the model
%   above?
