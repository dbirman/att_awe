function fit = fitCCBehavControlModel(adata,figs,model)
% CCBehavModel
%
% Fit the contrast (naka-rushton) and coherence (linear) models to the data
% obtained from the behavioral experiment. Use only the control condition
% data. This is just to compare linear vs. non-linear and constant vs.
% decreasing noise.

global fixedParams

adata = adata(~any(isnan(adata),2),:);
osize = size(adata,1);

fixedParams = struct;

adata = adata(adata(:,9)==-1,:);
disp(sprintf('Reducing data to %i control trials from %i',size(adata,1),osize));

%% Contrast Modeling Parameters
numParams = 0;
if strfind(model,'null')
    disp('(behavmodel) Fitting null contrast model');
    initparams.conslope = 0;
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
    initparams.conn = 1;
    numParams = numParams+3;
    initparams.conmodel = 2;
elseif strfind(model,'con-exp')
    disp('(behavmodel) Fitting exp contrast model');
    initparams.conalpha = [30 -inf inf];
    initparams.conkappa = [0.5 0 inf];
    initparams.conmodel = 3;
    numParams = numParams+2;
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
    numParams = numParams+3;
    initparams.cohn = 1;
    initparams.cohmodel = 2;
elseif strfind(model,'coh-exp')
    disp('(behavmodel) Fitting exp coherence model');
    initparams.cohmodel = 3;
    initparams.cohalpha = [30 -inf inf];
    initparams.cohkappa = [0.5 0 inf];
    numParams = numParams+2;
end

% freeze contrast and coherence at 1 so they force the other betas to
% similar values (i.e. sigma can't trade off with the other functions)
initparams.beta_control_con_conw = 1;
initparams.beta_control_con_cohw = [0 -inf inf];
initparams.beta_control_coh_cohw = 1;
initparams.beta_control_coh_conw = [0 -inf inf];
numParams = numParams+2;
    
if strfind(model,'poisson')
    disp('(behavmodel) Fitting poisson noise');
    initparams.poissonNoise = 1;
    initparams.sigma = 1;
else
    initparams.poissonNoise = 0;
    initparams.sigma = 1;
end

if strfind(model,'nobias')
    disp('(behavmodel) No bias');
    initparams.bias = 0;
else
    initparams.bias = [0 -inf inf];
    numParams = numParams+1;
end

if strfind(model,'stayswitch')
    disp('(behavmodel) Fitting four stay/switch parameters');
    initparams.right_correct = [0 -inf inf];
    initparams.right_incorr = [0 -inf inf];
    initparams.left_correct = [0 -inf inf];
    initparams.left_incorr = [0 -inf inf];
    numParams = numParams+4;
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

fit.modelstr = model;
fit.numParams = numParams;

function [bestparams,fit] = fitModel(params,adata,f,numParams)

[initparams, minparams, maxparams] = initParams(params);

options = optimoptions('fmincon','TolFun',0.05); % set a limit or it goes on foreeeeeeeeeeeever
bestparams = fmincon(@(p) fitBehavModel(p,adata,f),initparams,[],[],[],[],minparams,maxparams,[],options);

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

likelihood = 0;
% For each observation in adata, calculate log(likelihood) and sum

%     1       2         3        4      5      6     7      8       9
%   task - basecon - basecoh - conL - conR - cohL - cohR - resp - catch -
%      10      11        12
%   pedcon - pedcoh - correct

probs = zeros(size(adata,1),1);

for ai = 1:size(adata,1)
    obs = adata(ai,:);
    
    if ai>1
        prob = getObsProb(obs,params,adata(ai-1,:));
    else
        prob = getObsProb(obs,params,[]);
    end
    
    if prob==0
        warning('probably returned zero')
        prob = eps;
    end
    
    probs(ai) = prob;
    if prob >= 0
        likelihood = likelihood + log(prob);
    else
        warning('Probability returned non-useful value');
        keyboard
    end
end

likelihood = -likelihood;

% if likelihood==Inf
%     likelihood = 1000000;
% end

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
    title(sprintf('L: %2.3f.',likelihood));
end

function prob = getObsProb(obs,params,pobs)
%%
% check differences, adjust if necessary (to reflect actual visual amount shown)
if (obs(5)-obs(2))>(1-obs(2))
    keyboard;
end
if (obs(4)-obs(2))>(1-obs(2))
    keyboard;
end
if (obs(7)-obs(3))>(1-obs(3))
    keyboard;
end
if (obs(6)-obs(3))>(1-obs(3))
    keyboard;
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
            betas = [params.beta_att_coh_conw params.beta_att_coh_cohw];
        else
            % catch
            betas = [params.beta_unatt_coh_conw params.beta_unatt_coh_cohw];
        end
    case -2
        if obs(9)==0
            % main
            betas = [params.beta_att_con_conw params.beta_att_con_cohw];
        else
            betas = [params.beta_unatt_con_conw params.beta_unatt_con_cohw];
        end
end

if isfield(params,'right_correct') && ~isempty(pobs)
    extra = 0;
    if pobs(12)==1 && pobs(8)==1 % left correct
        extra = params.left_correct;
    elseif pobs(12)==1 && pobs(8)==2 % right correct
        extra = params.right_correct;
    elseif pobs(12)==0 && pobs(8)==1
        extra = params.left_incorr;
    elseif pobs(12)==0 && pobs(8)==2
        extra = params.right_incorr;
    end
    effect = betas * [conEff cohEff]' + params.bias + extra;
else
    effect = betas * [conEff cohEff]' + params.bias;
end

if params.poissonNoise
    if obs(8)==1
        prob = normcdf(0,effect,sqrt(abs(effect*params.sigma)),'upper');
    elseif obs(8)==0
        prob = normcdf(0,effect,sqrt(abs(effect*params.sigma)));
    else warning('failure'); keyboard
    end
else
    if obs(8)==1
        prob = normcdf(0,effect,params.sigma,'upper');
    elseif obs(8)==0
        prob = normcdf(0,effect,params.sigma);
    else warning('failure'); keyboard
    end
end

% if obs(8)==1, prob = 1-prob; end

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