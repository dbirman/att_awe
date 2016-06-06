function fit = fitCCBehavModel(adata,figs,model)
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

%% Contrast Modeling Parameters
numParams = 0;
if strfind(model,'null')
    disp('(behavmodel) Fitting null contrast model');
    initparams.conslope = 0;
    numParams = numParams+1;
    initparams.conmodel = 1;
elseif strfind(model,'con-linear')
    disp('(behavmodel) Fitting linear contrast model');
    initparams.conslope = [1 -inf inf];
    numParams = numParams+1;
    initparams.conmodel = 1;
elseif strfind(model,'con-naka')
    disp('(behavmodel) Fitting naka contrast model');
    initparams.conRmax = [30 -inf inf];
    initparams.conc50 = [0.75 0 1];
    initparams.conn = 2;
    numParams = numParams+2;
    initparams.conmodel = 2;
end
if strfind(model,'null')
    disp('(behavmodel) Fitting null coherence model');
    initparams.cohslope = 0;
    numParams = numParams+1;
    initparams.cohmodel = 1;
elseif strfind(model,'coh-linear')
    disp('(behavmodel) Fitting linear coherence model');
    initparams.cohslope = [10 -inf inf];
    numParams = numParams+1;
    initparams.cohmodel = 1;
elseif strfind(model,'coh-naka')
    disp('(behavmodel) Fitting naka coherence model');
    initparams.cohRmax = [1 -inf inf];
    initparams.cohc50 = [0.5 0 1];
    numParams = numParams+2;
    initparams.cohn = 1;
    initparams.cohmodel = 2;
end

initparams.scale = 1;

if strfind(model,'noalpha')
    initparams.alphacon = 1;
    initparams.alphacoh = 1;
    initparams.alphacon_att = 1;
    initparams.alphacoh_att = 1;
    initparams.alphacon_un = 1;
    initparams.alphacoh_un = 1;
else
    initparams.alphacon = [0.9 0 1];
    initparams.alphacoh = [0.9 0 1];
    initparams.alphacon_att = [0.8 0 1];
    initparams.alphacoh_att = [0.8 0 1];
    initparams.alphacon_un = [0.8 0 1];
    initparams.alphacoh_un = [0.8 0 1];
    numParams = numParams+6;
end

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
    initparams.sigma = 1.25;
else
    initparams.poissonNoise = 0;
    initparams.sigma = 1;
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
params = getParams(params);

params.conRmax = params.conRmax * params.scale;
params.cohslope = params.cohslope * params.scale;

% validate params
if params.alphacon < 0, params.alphacon = 0; end
if params.alphacon > 1, params.alphacon = 1; end
if params.alphacoh <0, params.alphacoh = 0; end
if params.alphacoh>1, params.alphacoh = 1; end

likelihood = 0;
% For each observation in adata, calculate log(likelihood) and sum
%   task - basecon - basecoh - conL - conR - cohL - cohR - resp - catch

probs = zeros(size(adata,1),1);

if f>0
    figure(f)
    clf
    hold on
%     title(sprintf('c50: %0.2f slope: %1.2f alphacon %0.2f alphacoh %0.2f',params.c50,params.slope,params.alphacon,params.alphacoh));
end

for ai = 1:size(adata,1)
    
    prob = -10;
    obs = adata(ai,:);
    
    prob = getObsProb(obs,params);
    
    probs(ai) = prob;
    if prob >= 0
        likelihood = likelihood + log(prob);
    else
        stop = 1;
    end
end

likelihood = -likelihood;

if f>0
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
    title(sprintf('Rmax: %2.2f c50: %2.2f Slope: %2.2f Sigma: %2.2f L: %2.3f.',params.conRmax,params.conc50,params.cohslope,params.sigma,likelihood));
end

function prob = getObsProb(obs,params)
%%

if obs(9)==1
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

    if params.unattNoise
        if obs(1)==-1
            params.sigma = params.sigma * params.conunatt;
        elseif obs(1)==-2
            params.sigma = params.sigma * params.cohunatt;
        end
    end
end

conEff = (conModel(obs(5),params)-conModel(obs(2),params)) - (conModel(obs(4),params)-conModel(obs(2),params));
cohEff = (cohModel(obs(7),params)-cohModel(obs(3),params)) - (cohModel(obs(6),params)-cohModel(obs(3),params));

if params.poissonNoise
    conProb = normcdf(0,conEff,sqrt(abs(conEff*params.sigma)));
    cohProb = normcdf(0,cohEff,sqrt(abs(cohEff*params.sigma)));
else
    conProb = normcdf(0,conEff,params.sigma);
    cohProb = normcdf(0,cohEff,params.sigma);
end

switch obs(1) % switch condition
    case 1
        % coherence control
        prob = cohProb * params.alphacoh + conProb * (1-params.alphacoh);
    case 2
        prob = conProb * params.alphacon + cohProb * (1-params.alphacon);
    case -1
        if obs(9)==0
            % main
            prob = cohProb * params.alphacoh_att + conProb *(1-params.alphacoh_att);
        else
            % catch
            prob = conProb * params.alphacon_un + cohProb * (1-params.alphacon_un);
        end
    case -2
        if obs(9)==0
            % main
            prob = conProb * params.alphacon_att + cohProb * (1-params.alphacon_att);
        else
            prob = cohProb * params.alphacoh_un + conProb * (1-params.alphacoh_un);
        end
end

if obs(8), prob = 1-prob; end
    
function out = conModel(con,params)
if params.conmodel==1
    out = params.conslope .* con;
elseif params.conmodel==2
    params.conn = round(params.conn);
    out = params.conRmax .* ((con.^params.conn) ./ (con.^params.conn + params.conc50.^params.conn));
end

function out = cohModel(coh,params)
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