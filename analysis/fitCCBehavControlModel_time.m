function fit = fitCCBehavControlModel_time(adata,figs,model)
% CCBehavModel
%
% Fit the contrast (naka-rushton) and coherence (linear) models to the data
% obtained from the behavioral experiment. Use only the control condition
% data. This is just to compare linear vs. non-linear and constant vs.
% decreasing noise.

global fixedParams

adata = adata(~any(isnan(adata),2),:);
osize = size(adata,1);

adata = adata(adata(:,9)==-1,:);
disp(sprintf('Reducing data to %i control trials from %i',size(adata,1),osize));

if any(adata(:,13)==1000)
    adata(:,13) = adata(:,13)/250;
end

%% Special condition: just getting BIC for a model
if isstruct(model)
    likelihood = fitBehavModel(model.params,adata,-1);
    fit.likelihood = likelihood;
    fit.BIC = 2*fit.likelihood + model.numParams * log(size(adata,1));
    return
elseif strfind(model,'sigma')
    initparams.conmodel = 4;
    initparams.cohmodel = 4;
    initparams.beta_control_con_conw = 1;
    initparams.beta_control_con_cohw = 0;
    initparams.beta_control_coh_cohw = 1;
    initparams.beta_control_coh_conw = 0;
    initparams.bias = 0;
    if strfind(model,'poisson')
        initparams.poissonNoise = 1;
        initparams.sigma = [0.002 eps inf];
    else
        initparams.poissonNoise = 0;
        initparams.sigma = [0.02 eps inf];
    end
    [~, fit] = fitModel(initparams,adata,-1,1);
    return
end

fixedParams = struct;
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

if strfind(model,'notime')
    initparams.time = 0;
else
    initparams.time = [0.5 -inf inf]; % this is the time gain which is the exponent on the time (1:2:4)
    numParams = numParams+1;
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

% 
if length(params.sigma)>1
    options = optimoptions('fmincon','TolFun',0.01); % set a limit or it goes on foreeeeeeeeeeeever
else
    options = optimoptions('fmincon'); % set a limit or it goes on foreeeeeeeeeeeever
end
bestparams = fmincon(@(p) fitBehavModel(p,adata,f),initparams,[],[],[],[],minparams,maxparams,[],options);

fit.params = getParams(bestparams);
[fit.likelihood] = fitBehavModel(bestparams,adata,0);
fit.BIC = 2*fit.likelihood + numParams * log(size(adata,1));

function likelihood = fitBehavModel(params,adata,f)
%%
if ~isstruct(params) && any(isnan(params))
    likelihood = Inf;
    return
end
if ~isstruct(params)
    params = getParams(params);
end

likelihood = 0;
% For each observation in adata, calculate log(likelihood) and sum

%     1       2         3        4      5      6     7      8       9
%   task - basecon - basecoh - conL - conR - cohL - cohR - resp - catch -
%      10      11        12
%   pedcon - pedcoh - correct

probs = zeros(size(adata,1),1);

% compute betas
betas = zeros(6,2);
betas(1,:) = [params.beta_control_coh_conw params.beta_control_coh_cohw];
betas(2,:) = [params.beta_control_con_conw params.beta_control_con_cohw];
% betas(3,:) = [params.beta_att_coh_conw params.beta_att_coh_cohw];
% betas(4,:) = [params.beta_unatt_coh_conw params.beta_unatt_coh_cohw];
% betas(5,:) = [params.beta_att_con_conw params.beta_att_con_cohw];
% betas(6,:) = [params.beta_unatt_con_conw params.beta_unatt_con_cohw];

% compute effects
conEffL = (conModel(adata(:,5),params)-conModel(adata(:,2),params)).*(adata(:,13).^params.time);;
conEffR = (conModel(adata(:,4),params)-conModel(adata(:,2),params)).*(adata(:,13).^params.time);;
conEff = conEffL - conEffR;
cohEffL = (cohModel(adata(:,7),params)-cohModel(adata(:,3),params)).*(adata(:,13).^params.time);;
cohEffR = (cohModel(adata(:,6),params)-cohModel(adata(:,3),params)).*(adata(:,13).^params.time);;
cohEff = cohEffL - cohEffR;

for ai = 1:size(adata,1)
    obs = adata(ai,:);
    
    if ai>1
        prob = getObsProb(obs,params,adata(ai-1,:),betas,conEff(ai),cohEff(ai),[conEffL(ai) conEffR(ai)],[cohEffL(ai) cohEffR(ai)]);
    else
        prob = getObsProb(obs,params,[],betas,conEff(ai),cohEff(ai),[conEffL(ai) conEffR(ai)],[cohEffL(ai) cohEffR(ai)]);
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
    plot(x,fcon*2^params.time,'--','Color',clist(1,:));
    plot(x,fcon*4^params.time,'-.','Color',clist(1,:));
    plot(x,fcoh,'Color',clist(3,:));
    plot(x,fcoh*2^params.time,'--','Color',clist(3,:));
    plot(x,fcoh*4^params.time,'-.','Color',clist(3,:));
    % now plot the unattended curves    
    title(sprintf('L: %2.3f.',likelihood));
end

function prob = getObsProb(obs,params,pobs,betas,conEff,cohEff,cons,cohs)
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
% MOVING OUTSIDE LOOP
% conEff = (conModel(obs(5),params)-conModel(obs(2),params)) - (conModel(obs(4),params)-conModel(obs(2),params));
% cohEff = (cohModel(obs(7),params)-cohModel(obs(3),params)) - (cohModel(obs(6),params)-cohModel(obs(3),params));

switch obs(1) % switch condition
    case 1
        % coherence control
        beta = betas(1,:);
    case 2
        % contrast control
        beta = betas(2,:);
    case -1
        if obs(9)==0
            % main coherence
            beta = betas(3,:);
        else
            % catch
            beta = betas(4,:);
        end
    case -2
        if obs(9)==0
            % main
            beta = betas(5,:);
        else
            beta = betas(6,:);
        end
end

extra = 0;
if isfield(params,'right_correct') && ~isempty(pobs)
    if pobs(12)==1 && pobs(8)==1 % left correct
        extra = params.left_correct;
    elseif pobs(12)==1 && pobs(8)==2 % right correct
        extra = params.right_correct;
    elseif pobs(12)==0 && pobs(8)==1
        extra = params.left_incorr;
    elseif pobs(12)==0 && pobs(8)==2
        extra = params.right_incorr;
    end
end
effect = beta * [conEff cohEff]' + params.bias + extra;

if params.poissonNoise
    noise = sqrt(abs(sum([mean(cons) mean(cohs)])));
    if obs(8)==1
        prob = normcdf(0,effect,noise*params.sigma,'upper');
    elseif obs(8)==0
        prob = normcdf(0,effect,noise*params.sigma);
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