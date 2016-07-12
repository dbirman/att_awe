function fit = fitCCBehavModelROI(adata,figs,model,subj)
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

fixedParams = struct;

%% Load subj data
if strfind(getenv('OS'),'Windows')
    load(fullfile(sprintf('C:/Users/Dan/proj/COHCON_DATA/%s_fitroi.mat',subj)));
else
    load(sprintf('~/data/cohcon_localizer/%s_fitroi.mat',subj));
end
initparams.conRmax = mean([fitroi.roiparams{1}.Rmax fitroi.roiparams{2}.Rmax]);
initparams.conc50 = mean([fitroi.roiparams{1}.c50 fitroi.roiparams{2}.c50]);
initparams.conn = 1;
initparams.conmodel = 2;
initparams.cohslope = mean([fitroi.roiparams{19}.slope fitroi.roiparams{20}.slope]);
initparams.cohmodel = 1;

%% Setup params
numParams = 0;

initparams.scale = 1;

% beta parameters
initparams.beta_control_con_conw = 1;
initparams.beta_control_con_cohw = [0 -inf inf];
initparams.beta_control_coh_conw = [0 -inf inf];
initparams.beta_control_coh_cohw = 1;
initparams.beta_att_con_conw = 1;
initparams.beta_att_con_cohw = [0 -inf inf];
initparams.beta_att_coh_conw = [0 -inf inf];
initparams.beta_att_coh_cohw = 1;
initparams.beta_unatt_con_conw = [0 -inf inf];
initparams.beta_unatt_con_cohw = 1;
initparams.beta_unatt_coh_conw = 1;
initparams.beta_unatt_coh_cohw = [0 -inf inf];

numParams = numParams+6;

if strfind(model,'nounatt')
    initparams.conunatt = 1;
    initparams.cohunatt = 1;
    initparams.unattNoise = 0;
elseif strfind(model,'unattnoise')
    initparams.conunatt = [0.9 0 1];
    initparams.cohunatt = [0.9 0 1];
    numParams = numParams+2;
    initparams.unattNoise = 1;
else
    initparams.conunatt = [0.9 0 1];
    initparams.cohunatt = [0.9 0 1];
    numParams = numParams+2;
    initparams.unattNoise = 0;
end

if strfind(model,'poisson')
    initparams.poissonNoise = 1;
    initparams.sigma = [1.25 0 inf];
else
    initparams.poissonNoise = 0;
    initparams.sigma = [1 0 inf];
end
% initparams.sigma = 0.01088;
% initparams.poissonNoise = 0;
    
if strfind(model,'nobias')
    initparams.bias = 0;
else
    initparams.bias = [0 -inf inf];
    numParams = numParams+1;
end

if strfind(model,'onlycatch')
    fixedParams.onlycatch = 1;
elseif strfind(model,'nocatch')
    fixedParams.onlycatch = 0;
else
    fixedParams.onlycatch = -1;
end

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
if params.conunatt<0, params.conunatt = 0; end
if params.cohunatt<0, params.cohunatt = 0;  end

if params.conmodel==2
    params.conRmax = params.conRmax * params.scale;
else
    params.conslope = params.conslope * params.scale;
end
if params.cohmodel==2
    params.cohRmax = params.cohRmax * params.scale;
else
    params.cohslope = params.cohslope * params.scale;
end

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
global fixedParams

if obs(9)==1 && fixedParams.onlycatch==0
    % catch trial, and we don't want any catch trials
    prob = 1; return
elseif obs(9)==-1 && fixedParams.onlycatch==1
    % not a catch trial, and we want only catch trials
    prob = 1; return
end

if obs(9)==1
    if params.unattNoise
        if obs(1)==-1
            params.sigma = params.sigma * params.conunatt;
        elseif obs(1)==-2
            params.sigma = params.sigma * params.cohunatt;
        end
    else     
        if obs(1)==-1
            % THIS IS A CATCH TRIAL IN A COHERENCE RUN, ADJUST CONTRAST
            if params.conmodel==1
                params.conslope = params.conslope*params.conunatt;
            else
                params.conRmax = params.conRmax*params.conunatt;
            end
        elseif obs(1)==-2
            % ADJUST COHERENCE
            if params.cohmodel==1
                params.cohslope = params.cohslope * params.cohunatt;
            else
                params.cohRmax = params.cohRmax*params.cohunatt;
            end
        end
    end
end

conEff = (conModel(obs(5),params)-conModel(obs(2),params)) - (conModel(obs(4),params)-conModel(obs(2),params));
cohEff = (cohModel(obs(7),params)-cohModel(obs(3),params)) - (cohModel(obs(6),params)-cohModel(obs(3),params));

switch obs(1) % switch condition
    case 1
        % coherence control
        betas = [params.beta_control_coh_conw params.beta_control_coh_cohw];
    case 2
        % contrast control
        betas = [params.beta_control_con_conw params.beta_control_con_cohw];
    case -1
        if obs(9)==0
            % main coherence
            betas = [params.beta_control_coh_conw params.beta_att_coh_cohw];
        else
            % catch
            betas = [params.beta_control_coh_conw params.beta_unatt_coh_cohw];
        end
    case -2
        if obs(9)==0
            % main
            betas = [params.beta_control_con_conw params.beta_att_con_cohw];
        else
            betas = [params.beta_control_con_conw params.beta_unatt_con_cohw];
        end
end
effect = betas * [conEff cohEff]' + params.bias;

if params.poissonNoise
    prob = normcdf(0,effect,sqrt(abs(effect*params.sigma)));
else
    prob = normcdf(0,effect,params.sigma);
end

if obs(8), prob = 1-prob; end
    
function out = conModel(con,params)
if con==0, out=0; return; end
if params.conmodel==1
    out = params.conslope .* con;
elseif params.conmodel==2
    params.conn = round(params.conn);
    out = params.conRmax .* ((con.^params.conn) ./ (con.^params.conn + params.conc50.^params.conn));
end

function out = cohModel(coh,params)
if coh==0, out=0; return; end
if params.cohmodel==1
    out = params.cohslope .* coh;
elseif params.cohmodel==2
    params.cohn = round(params.cohn);
    out = params.cohRmax .* ((coh.^params.cohn) ./ (coh.^params.cohn + params.cohc50.^params.cohn));
end

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