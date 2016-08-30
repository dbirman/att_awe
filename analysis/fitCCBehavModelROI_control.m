function fit = fitCCBehavModelROI_control(adata,figs,model,subj,subjb)
% CCBehavModel
%
% Fit the contrast (naka-rushton) and coherence (linear) models to the data
% obtained from the behavioral experiment. Roughly we will do the
% following:
%
% Model Types:
% 'null' : model with no effects Rmax=0 slope = 0;
% 'con-linear'
% 'coh-linear'
% 'con-naka'
% 'coh-naka'
% 'con-n'
% 'coh-n'
global fixedParams

adata = adata(~any(isnan(adata),2),:);
osize = size(adata,1);
adata = adata(adata(:,9)==-1,:);
disp(sprintf('Reducing data to %i control trials from %i',size(adata,1),osize));

fixedParams = struct;

%% Load subj data
if strfind(getenv('OS'),'Windows')
    load(fullfile(datafolder,sprintf('%s_fitroi.mat',subj)));
else
    load(sprintf('~/data/cohcon_localizer/%s_fitroi.mat',subj));
end

x = 0:.001:2;
con_lv1 = conModel(x,fitroi.roiparams{1})-(fitroi.roiparams{1}.offset-fitroi.roiparams{1}.conalpha);
con_rv1 = conModel(x,fitroi.roiparams{2})-(fitroi.roiparams{2}.offset-fitroi.roiparams{2}.conalpha);
coh_lmt = cohModel(x,fitroi.roiparams{19})-(fitroi.roiparams{19}.offset-fitroi.roiparams{19}.cohalpha);
coh_rmt = cohModel(x,fitroi.roiparams{2})-(fitroi.roiparams{20}.offset-fitroi.roiparams{20}.cohalpha);

fixedParams.x = x;
fixedParams.con = mean([con_lv1;con_rv1]);
fixedParams.coh = mean([coh_lmt;coh_rmt]);

%% Setup params
numParams = 0;

initparams.conmodel=4;
initparams.cohmodel=4;

initparams.scale = 1;

% beta parameters

% we are going to force the model to use the same values that come from the
% behavior and only fit sigma, this is better than allowing it to vary
fields = {'beta_control_con_conw','beta_control_con_cohw','beta_control_coh_conw','beta_control_coh_cohw','bias'};
load(fullfile(datafolder,sprintf('%s_data.mat',subjb)));
for fi = 1:length(fields)
    initparams.(fields{fi}) = data.fits{1}.params.(fields{fi});
end
% initparams.beta_control_con_conw = 1;
% initparams.beta_control_con_cohw = [0 -inf inf];
% initparams.beta_control_coh_conw = [0 -inf inf];
% initparams.beta_control_coh_cohw = 1;

% numParams = numParams+2;

initparams.sigma = [1/50 0 inf];
initparams.poissonNoise = 0;

numParams = numParams+1;
    
% if strfind(model,'nobias')
%     initparams.bias = 0;
% else
%     initparams.bias = [0 -inf inf];
%     numParams = numParams+1;
% end

%% Prep and Call
if ieNotDefined('figs')
    figs = 0;
end

% Call fmins
if figs
    f = figure;
else
    f = -1;
end

[~, fit] = fitModel(initparams,adata,f,numParams);

function [bestparams,fit] = fitModel(params,adata,f,numParams)

[initparams, minparams, maxparams] = initParams(params);

bestparams = fmincon(@(p) fitBehavModel(p,adata,f),initparams,[],[],[],[],minparams,maxparams);

fit.params = getParams(bestparams);
[fit.likelihood] = fitBehavModel(bestparams,adata,0);
fit.BIC = 2*fit.likelihood + numParams * log(size(adata,1));

function likelihood = fitBehavModel(params,adata,f)
%%
if any(isnan(params))
    likelihood = Inf;
    return
end
params = getParams(params);

if params.sigma<0.0000001, params.sigma=0.0000001; end

likelihood = 0;
% For each observation in adata, calculate log(likelihood) and sum
%   task - basecon - basecoh - conL - conR - cohL - cohR - resp - catch

probs = zeros(size(adata,1),1);

for ai = 1:size(adata,1)
    obs = adata(ai,:);
    
    prob = getObsProb(obs,params);
    
    probs(ai) = prob;
    if prob >= 0
        likelihood = likelihood + log(prob);
    else
        warning('Probability returned non-useful value');
        keyboard
    end
end

likelihood = -likelihood;

if likelihood<.001
    warning('Potential failure...');
    keyboard
end

if f>0
    warning('No figure implemented');
    keyboard
    figure(f)
    clf
    hold on
    clist = brewermap(3,'PuOr');
    x = 0:.01:1;
    fcon = conModel(x,params);
%     fconp = 1-normcdf(0,fcon,params.sigma_con);
    fcoh = cohModel(x,params);
%     fcohp = 1-normcdf(0,fcoh,params.sigma_coh);
    plot(x,fcon,'Color',clist(1,:));
    plot(x,fcoh,'Color',clist(3,:));
    % now plot the unattended curves    
    if params.conmodel==1
        params.conslope = params.conslope*params.conunatt;
    else
        params.conRmax = params.conRmax*params.conunatt;
    end
    if params.cohmodel==1
        params.cohslope = params.cohslope * params.cohunatt;
    else
        params.cohRmax = params.cohRmax*params.cohunatt;
    end
    fcon = conModel(x,params);
    fcoh = cohModel(x,params);
    plot(x,fcon,'--','Color',clist(1,:));
    plot(x,fcoh,'--','Color',clist(3,:));
    title(sprintf('L: %2.3f.',likelihood));
end

function prob = getObsProb(obs,params)
%%
% global fixedParams

% basecon = fixedParams.con(find(fixedParams.x>=obs(2),1));
% lcon = fixedParams.con(find(fixedParams.x>=obs(5),1));
% rcon = fixedParams.con(find(fixedParams.x>=obs(4),1));
% conEff = (lcon-basecon) - (rcon-basecon);
% basecoh = fixedParams.coh(find(fixedParams.x>=obs(3),1));
% lcoh = fixedParams.coh(find(fixedParams.x>=obs(7),1));
% rcoh = fixedParams.coh(find(fixedParams.x>=obs(6),1));
% cohEff = (lcoh-basecoh) - (rcoh-basecoh);
conEff = (conModel(obs(5),params)-conModel(obs(2),params)) - (conModel(obs(4),params)-conModel(obs(2),params));
cohEff = (cohModel(obs(7),params)-cohModel(obs(3),params)) - (cohModel(obs(6),params)-cohModel(obs(3),params));

switch obs(1) % switch condition
    case 1
        % coherence control
        betas = [params.beta_control_coh_conw params.beta_control_coh_cohw];
    case 2
        % contrast control
        betas = [params.beta_control_con_conw params.beta_control_con_cohw];
end
effect = betas * [conEff cohEff]' + params.bias*params.sigma;

prob = normcdf(0,effect,params.sigma);

if obs(8), prob = 1-prob; end

function [initparams, minparams, maxparams] = initParams(params)

global fixedParams

fixedParams.strs = fields(params);
fixedParams.num = length(fixedParams.strs);

initparams = [];
minparams = [];
maxparams = [];
indexes = zeros(1,fixedParams.num);
count = 1;

fixed = zeros(1,fixedParams.num);
optim = zeros(1,fixedParams.num);

for i = 1:fixedParams.num
    cvals = params.(fixedParams.strs{i});
    
    if length(cvals)==1
        fixedParams.(fixedParams.strs{i}) = cvals;
        fixed(i) = 1;
    elseif length(cvals)==3
        initparams = [initparams cvals(1)];
        minparams = [minparams cvals(2)];
        maxparams = [maxparams cvals(3)];
        indexes(i) = count;
        count = count+1;
    elseif length(cvals)==2 || length(cvals)>3
        % optimizer
        fixedParams.(fixedParams.strs{i}) = cvals;
        optim(i) = 1;
    else
        error('You initialized a parameter with the wrong initial values... unable to interpret');
    end
end
fixedParams.optim = optim;
fixedParams.fixed = fixed;
fixedParams.idx = indexes;

function p = getParams(params)

global fixedParams

for i = 1:fixedParams.num
    if fixedParams.fixed(i)
        p.(fixedParams.strs{i}) = fixedParams.(fixedParams.strs{i});
    else
        p.(fixedParams.strs{i}) = params(fixedParams.idx(i));
    end
end