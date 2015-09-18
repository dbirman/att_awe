function simData = cc_runBehavSimulations()

%% This is a simulation function to compare different possible models.

%% Model A

% Model A is the following structure:
% I(L,con) = Stimulus Contrast Intensity on Left Side
% S(L,con) = f(I) naka-rushton eqn
% D(R,con) = S(R,con)-S(L,con) diff function
% P(R) = f(a * D(con) + (1-a) * D(coh)) | a is attention to contrast, f is
% the weibull function

s = struct;

%% Naka-Ruhston

tasks = {'coherence','contrast'};
params = {'N','p','q','sigma'};
values = {{3,1,1,0.25},{2,1,1,0.25}};

for ti = 1:length(tasks)
    for pi = 1:length(params)
        s.(tasks{ti}).(params{pi}) = values{ti}{pi};
    end
end

% Setup S = f(I)

for ti = 1:length(tasks)
    s.(tasks{ti}).naka = cc_naka(0:.01:1,s.(tasks{ti}));
end

% Plot Naka

s.map = brewermap(15,'PuOr');
s.contrast.color = s.map(5,:); s.coherence.color = s.map(11,:);

figure
hold on
legs = {};
for ti = 1:length(tasks)
    plot(0:.01:1,s.(tasks{ti}).naka,'-','Color',s.(tasks{ti}).color);
    legs{end+1} = tasks{ti};
end
legend(legs);

%% Generate Randomized Trials

trials = 5000;
tasks = [1 2];
% calculate basic values
conSide = datasample(tasks,trials);
cohSide = datasample(tasks,trials);
task = datasample(tasks,trials);
side = zeros(1,trials);
rCon = zeros(1,trials);
lCon = zeros(1,trials);
rCoh = zeros(1,trials);
lCoh = zeros(1,trials);
for i = 1:trials
    if task(i)==1
        side(i) = cohSide(i);
    else
        side(i) = conSide(i);
    end
    if conSide(i)==1
        % con left
        rCon(i) = .6;
        lCon(i) = .6 + rand*.35;
    else
        lCon(i) = .6;
        rCon(i) = .6 + rand*.35;
    end
    if cohSide(i) == 1
        rCoh(i) = .1;
        lCoh(i) = .1 + rand*.85;
    else
        lCoh(i) = .1;
        rCoh(i) = .1 + rand*.85;
    end
end

%% Calculate Strength for Trials

% calculate S values
lConS = interpNaka(lCon,s,'contrast');
rConS = interpNaka(rCon,s,'contrast');
lCohS = interpNaka(lCoh,s,'coherence');
rCohS = interpNaka(rCoh,s,'coherence');

%% Calculate Diff

dSCon = rConS - lConS;
dSCoh = rCohS - lCohS;
dCon = rCon-lCon;
dCoh = rCoh-lCoh;

%% Calculate Choice
con = .025; % effect of coherence on contrast
coh = .4; % effect of contrast on coherence

side = zeros(1,trials);
for ti = 1:trials
    if task(ti)==1
        % coherence main
        side(ti) = (dSCoh(ti)*(1-coh)+dSCon(ti)*coh) > 0;
    else
        side(ti) = (dSCoh(ti)*(con)+dSCon(ti)*(1-con)) > 0;
    end
end
%% Plotz

conbins = -.3:.05:.3; conrange = -.275:.05:.275;
cohbins = -.8:.1:.8; cohrange = -.75:.1:.75;
taskS = {'coherence','contrast'};
figure
for ti = 1:length(tasks)
    subplot(2,1,ti), hold on
    % one more beta value here, we have to deal with that
    % get actual data points
    conBinned = binData(side(task==ti),dCon(task==ti),conbins);
    conMean = cellfun(@mean,conBinned);
    conStd = cellfun(@(A) std(A)/sqrt(length(A))*1.96,conBinned);
    % coherence data points
    cohBinned = binData(side(task==ti),dCoh(task==ti),cohbins);
    cohMean = cellfun(@mean,cohBinned);
    cohStd = cellfun(@(A) std(A)/sqrt(length(A))*1.96,cohBinned);

    errorbar(conrange,conMean,conStd,'-','Color',s.contrast.color);
    errorbar(cohrange,cohMean,cohStd,'-','Color',s.coherence.color);
            
end
%% Regression

taskS = {'coherence','contrast'};
for ti = 1:length(tasks)
    subplot(2,1,ti), hold on
    X = [dCoh(task==ti)',dCon(task==ti)'];
    B = glmfit(X,side(task==ti)','binomial','link','logit');
    
    y = glmval(B,[zeros(61,1),(-.3:.01:.3)'],'logit');
    plot(-.3:.01:.3,y,'-','Color',s.contrast.color);
    y = glmval(B,[(-.85:.01:.85)',zeros(171,1)],'logit');
    plot(-.85:.01:.85,y,'-','Color',s.coherence.color)
    title(taskS{ti});
end

simData = s;

%% Behavioral Data Match

fullData = loadFullData;

%% Functions

% def: function w = weibull(x,params) params = T, Beta, lambda, min

function S = cc_interpNaka(I,s,task)

if length(I)>1
    for i = 1:length(I)
        S(i) = interpNaka(I(i),s,task);
    end
else
    S = interp1(s.(task).I,s.(task).naka,I) + s.(task).noise*randn;
end


