function fit = fitCCBehavControlModel_fmri(adata,figs,model,confit,cohfit,lapserate)
% CCBehavModel
%
% Copy of fitCCBehavControlModel that allows you to fit a doublesigma (one
% sigma for contrast, one for coherence)

global fixedParams

adata = adata(~any(isnan(adata),2),:);
osize = size(adata,1);

adata = adata(adata(:,9)==-1,:);
disp(sprintf('Reducing data to %i control trials from %i',size(adata,1),osize));

%% Lower basecontrast

%     1       2         3        4      5      6     7      8       9
%   task - basecon - basecoh - conL - conR - cohL - cohR - resp - catch -
%      10      11        12
%   pedcon - pedcoh - correct

if any(adata(:,2)>0)
    disp('BASE CONTRAST LOWERED TO ZERO');
    adata(:,4) = adata(:,4)-adata(:,2);
    adata(:,5) = adata(:,5)-adata(:,2);
    adata(:,2) = adata(:,2)-adata(:,2);
end

%% Special condition: just getting BIC for a model
if isstruct(model)
    likelihood = fitBehavModel(model.params,adata,-1);
    fit.likelihood = likelihood;
    fit.BIC = 2*fit.likelihood + model.numParams * log(size(adata,1));
    return
elseif strfind(model,'sigma')
    disp('Fitting sigma...');
    % SPECIAL CONDITION: Fitting sigma parameter
    fixedParams.x = 0:.001:1;
    if isstruct(confit)
        fixedParams.con = conModel(fixedParams.x,confit.params);
        fixedParams.coh = cohModel(fixedParams.x,cohfit.params);
    else
        fixedParams.con = confit;
        fixedParams.coh = cohfit;
    end
    initparams.conmodel = 4;
    initparams.cohmodel = 4;
    initparams.beta_control_con_conw = 1;
    initparams.beta_control_con_cohw = [0 -1 1];
    initparams.beta_control_coh_cohw = 1;
    initparams.beta_control_coh_conw = [0 -1 1];
    initparams.bias = [0 -1 1];
    if strfind(model,'poisson')
        initparams.poissonNoise = 1;
        if strfind(model,'doublesigma')
            initparams.sigmacon = [0.002 eps 1];
            initparams.sigmacoh = [0.002 eps 1];
        else
            initparams.sigma = [0.002 eps 1];
        end
    else
        initparams.poissonNoise = 0;
        if strfind(model,'doublesigma')
            initparams.sigmacon = [0.02 eps 1];
            initparams.sigmacoh = [0.02 eps 1];
        else
            initparams.sigma = [0.02 eps 1];
        end
    end
    initparams.lapse = lapserate;
    [~, fit] = fitModel(initparams,adata,-1);
    return
end

function [bestparams,fit] = fitModel(params,adata,f)

[initparams, minparams, maxparams] = initParams(params);

% 
if (isfield(params,'sigmacon') && length(params.sigmacon)>1) || (isfield(params,'sigma') && length(params.sigma)>1)
    options = optimoptions('fmincon','TolFun',0.01); % set a limit or it goes on foreeeeeeeeeeeever
else
    options = optimoptions('fmincon'); % set a limit or it goes on foreeeeeeeeeeeever
end
bestparams = fmincon(@(p) fitBehavModel(p,adata,f),initparams,[],[],[],[],minparams,maxparams,[],options);

fit.params = getParams(bestparams);
[fit.likelihood] = fitBehavModel(bestparams,adata,0);
fit.BIC = 2*fit.likelihood + length(bestparams) * log(size(adata,1));
fit.AIC = 2*fit.likelihood + length(bestparams) * 2;
fit.numParams = length(bestparams);

function likelihood = fitBehavModel(params,adata,f)
%%
if ~isstruct(params) && any(isnan(params))
    likelihood = Inf;
    return
end
if ~isstruct(params)
    params = getParams(params);
end

if isfield(params,'sigmacon') && params.sigmacon < eps
    likelihood = inf;
    return
end
if isfield(params,'sigma') && params.sigma < eps
    likelihood = inf;
    return
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
conEffL = (conModel(adata(:,4),params)-conModel(adata(:,2),params));
conEffR = (conModel(adata(:,5),params)-conModel(adata(:,2),params));
conEff = conEffR - conEffL;
cohEffL = (cohModel(adata(:,6),params)-cohModel(adata(:,3),params));
cohEffR = (cohModel(adata(:,7),params)-cohModel(adata(:,3),params));
cohEff = cohEffR - cohEffL;

for ai = 1:size(adata,1)
    obs = adata(ai,:);
    
    if ai>1
        prob = getObsProb(obs,params,adata(ai-1,:),betas,conEff(ai),cohEff(ai),[conEffL(ai) conEffR(ai)],[cohEffL(ai) cohEffR(ai)]);
    else
        prob = getObsProb(obs,params,[],betas,conEff(ai),cohEff(ai),[conEffL(ai) conEffR(ai)],[cohEffL(ai) cohEffR(ai)]);
    end
    
    % add lapse rate
    prob = params.lapse + (1-2*params.lapse)*prob;
    
    if prob==0
%         warning('probability returned zero')
        prob = eps;
    end
    
    probs(ai) = prob;
    if prob >= 0
        likelihood = likelihood + log(prob);
    elseif isnan(prob)
        warning('Probability returned non-useful value');
        prob = eps;
        likelihood = likelihood + log(prob);
%         keyboard
    end
end

likelihood = -likelihood;

if likelihood<.001
    warning('Potential failure...');
    keyboard
end

if 1
    figure(1)
    clf
    hold on
    clist = brewermap(3,'PuOr');
    x = 0:.001:1;
    if isfield(params,'sigmacon')
        fcon = params.sigmacon*conModel(x,params);
        fcoh = params.sigmacoh*cohModel(x,params);
    else
        fcon = params.sigma*conModel(x,params);
        fcoh = params.sigma*cohModel(x,params);
    end
    plot(x,fcon,'Color',clist(1,:));
    plot(x,fcoh,'Color',clist(3,:));
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

if isfield(params,'sigmacon')
    if obs(1)==1
        usesigma = params.sigmacoh;
    else
        usesigma = params.sigmacon;
    end
else
    usesigma = params.sigma;
end

if params.poissonNoise
    noise = sqrt(abs(sum(beta*[mean(cons) mean(cohs)]')));

    if obs(8)==1
        prob = normcdf(0,effect,noise*usesigma,'upper');
    elseif obs(8)==0
        prob = normcdf(0,effect,noise*usesigma);
    else warning('failure'); keyboard
    end
else
    if obs(8)==1
        prob = normcdf(0,effect,usesigma,'upper');
    elseif obs(8)==0
        prob = normcdf(0,effect,usesigma);
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