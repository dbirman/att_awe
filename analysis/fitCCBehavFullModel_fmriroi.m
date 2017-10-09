function fit = fitCCBehavFullModel_fmriroi(adata,sigmaval,model,confit,cohfit,lapserate,crossval)
% CCBehavModel
%
% Copy of fitCCBehavControlModel_fmri that allows you to fit an ROI variant
% instead of con/coh response functions, this has a contrast/coherence
% response for each ROI (no interactions) which are then weighted. ROI
% dropout can be specified. 

global fixedParams

adata = adata(~any(isnan(adata),2),:);
osize = size(adata,1);

%% Lower basecontrast

%     1       2         3        4      5      6     7      8       9
%   task - basecon - basecoh - conL - conR - cohL - cohR - resp - catch -
%      10      11        12
%   pedcon - pedcoh - correct

if any(adata(:,2)>0)
%     disp('BASE CONTRAST LOWERED TO ZERO');
    adata(:,4) = adata(:,4)-adata(:,2);
    adata(:,5) = adata(:,5)-adata(:,2);
    adata(:,2) = adata(:,2)-adata(:,2);
end

if crossval
%     disp('(fitCCBehavControlModel_fmri) CROSSVAL INITIATED');
    
    datapoints = 1:size(adata,1);
    num = floor(size(adata,1)/10);
    
    aprobs = []; trainr2 = []; testr2 = []; sigmas = []; aresp = []; pseudor2 = zeros(1,10); like = zeros(1,10); nulllike = zeros(1,10);
    for fold = 1:10
        if fold < 10
            test = datapoints((fold-1)*num+1:(fold-1)*num+num);
        else
            test = datapoints((fold-1)*num+1:end);
        end
        train = setdiff(datapoints,test);
        
        traindata = adata(train,:);
        testdata = adata(test,:);
        
        if ~isempty(intersect(train,test))
            warning('Failure to fully separate training and testing data');
            keyboard
        end
        
        trainfit = fitCCBehavControlModel_fmri(traindata,sigmaval,model,confit,cohfit,lapserate,0);
%         trainr2(fold) = trainfit.pseudor2;
        sigmas(fold) = trainfit.params.sigma;
        testfit = fitCCBehavControlModel_fmri(testdata,sigmaval,trainfit,confit,cohfit,lapserate,0);
%         testr2(fold) = testfit.pseudor2;
        aresp = [aresp ;testdata(:,8)];
        aprobs = [aprobs testfit.probs];
        pseudor2(fold) = testfit.pseudor2;
        like(fold) = testfit.likelihood;
        nulllike(fold) = testfit.null.likelihood;
    end
    
    fit = fitCCBehavControlModel_fmri(adata,sigmaval,model,confit,cohfit,lapserate,0);
    fit.cv.sigmas = sigmas;
    fit.cv.like = like;
    fit.cv.nulllike = nulllike;
    fit.cv.pseudor2 = 1-(sum(like)/sum(nulllike));
    fit.cv.resp = aresp;
    fit.cv.trainpr2 = trainr2;
    fit.cv.testpr2 = testr2;
    fit.cv.aprobs = aprobs;
    fit.cv.cd = nanmean(aprobs(aresp==1))-nanmean(1-aprobs(aresp==0));
    
    return
end

%% Special condition: just getting BIC for a model
if isstruct(model)
    fit = struct;
    fit.null = fitCCBehavControlModel_fmri(adata,sigmaval,'null',confit,cohfit,0,0);
    [fit.likelihood,probs] = fitBehavModel(model.params,adata,-1);
    fit.BIC = 2*fit.likelihood + model.numParams * log(size(adata,1));
    fit.AIC = 2*fit.likelihood + model.numParams * 2;
    fit.params = model.params;
    fit.probs = probs;
    fit.resp = adata(:,8);
    fit.cd = nanmean(probs(fit.resp==1))-nanmean(1-probs(fit.resp==0));
    fit.pseudor2 = 1 - (fit.likelihood / fit.null.likelihood);
    return
elseif strfind(model,'null')
    % ONLY ALLOW BIAS TO CHANGE (intercept only model)
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
    initparams.beta_control_con_conw = 0;
    initparams.beta_control_con_cohw = 0;%[0 -1 1];
    initparams.beta_control_coh_cohw = 0;
    initparams.beta_control_coh_conw = 0;%[0 -1 1];
    initparams.bias = [0 -1 1];
    initparams.sigma = 1;
    initparams.coh_gain = 0;
    initparams.poissonNoise = 0;
    initparams.lapse = 0;
    [~, fit] = fitModel(initparams,adata,-1);
    return
elseif strfind(model,'freeze')
    % ONLY ALLOW SIGMA TO CHANGE
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
    initparams.beta_control_con_cohw = 0;%[0 -1 1];
    initparams.beta_control_coh_cohw = 1;
    initparams.beta_control_coh_conw = 0;%[0 -1 1];
    initparams.bias = 0;%;[0 -1 1];
    initparams.sigma = [0.1 eps 1];
    initparams.poissonNoise = 0;
    initparams.lapse = lapserate;
    [~, fit] = fitModel(initparams,adata,-1);
    return
elseif strfind(model,'sigma')
%     disp('Fitting sigma...');
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
    if strfind(model,'gain')
        initparams.coh_gain = [1.2 0 inf];
        initparams.sigma = sigmaval.sigma; % we're doing something funky here by using figs, but fuckit
        initparams.poissonNoise = 0;
        initparams.beta_control_con_cohw = sigmaval.beta_control_con_cohw;
        initparams.beta_control_coh_conw = sigmaval.beta_control_coh_conw;
        initparams.bias = sigmaval.bias;
    else
        initparams.coh_gain = 1;
        if strfind(model,'poisson')
            initparams.poissonNoise = 1;
            if strfind(model,'doublesigma')
%                 initparams.sigmacon = [0.1 eps 1];
%                 initparams.sigmacoh = [0.1 eps 1];
            else
                initparams.sigma = [sigmaval eps 1];
            end
        elseif strfind(model,'mergenoise')
            % mixture of poisson and non-poisson noise
            initparams.psigma = [.17 eps 1];
            initparams.sigma = [.11 eps 1];
            initparams.merge = [0.5 0 1];
        else
            initparams.poissonNoise = 0;
            if strfind(model,'doublesigma')
%                 initparams.sigmacon = [0.1 eps 1];
%                 initparams.sigmacoh = [0.1 eps 1];
            else
                initparams.sigma = [sigmaval eps 1];
            end
        end
    end
    
    if strfind(model,'stayswitch')
%         disp('(behavmodel) Fitting four stay/switch parameters');
        initparams.right_correct = [0 -inf inf];
        initparams.right_incorr = [0 -inf inf];
        initparams.left_correct = [0 -inf inf];
        initparams.left_incorr = [0 -inf inf];
    end
    initparams.lapse = lapserate;
    [~, fit] = fitModel(initparams,adata,-1);
    return
end

function [bestparams,fit] = fitModel(params,adata,f)

[initparams, minparams, maxparams] = initParams(params);

options = optimoptions('fmincon','Algorithm','active-set','TolFun',1,'TolCon',1,'Display','off'); % set a limit or it goes on foreeeeeeeeeeeever

bestparams = fmincon(@(p) fitBehavModel(p,adata,f),initparams,[],[],[],[],minparams,maxparams,[],options);

fit.params = getParams(bestparams);
[fit.likelihood, probs] = fitBehavModel(bestparams,adata,0);
fit.BIC = 2*fit.likelihood + length(bestparams) * log(size(adata,1));
fit.AIC = 2*fit.likelihood + length(bestparams) * 2;
fit.probs = probs;
fit.resp = adata(:,8);
fit.cd = nanmean(probs(fit.resp==1))-nanmean(1-probs(fit.resp==0));
fit.numParams = length(bestparams);

function [likelihood, probs] = fitBehavModel(params,adata,f)
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
if params.coh_gain ~= 1
    cohEffL = cohEffL * params.coh_gain;
    cohEffR = cohEffR * params.coh_gain;
end
cohEff = cohEffR - cohEffL;

probs = zeros(1,size(adata,1));

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
    
    if isnan(prob)
        warning('Probability returned non-useful value');
        prob = eps;
        keyboard
    end
    probs(ai) = prob;
end

likelihood = -sum(log(probs));

if likelihood<.001
    warning('Potential failure...');
    keyboard
end

if 0
    figure(1)
    clf
    resp = adata(:,8);
    subplot(211)
    bins = 0.05:.1:.95;
    hist(probs(resp==1),bins);
    title('Response RIGHT');
    subplot(212);
    hist(1-probs(resp==0),bins);
    title('Response LEFT');
    stop = 1;
%     clf
%     hold on
%     clist = brewermap(3,'PuOr');
%     x = 0:.001:1;
%     if isfield(params,'sigmacon')
%         fcon = params.sigmacon*conModel(x,params);
%         fcoh = params.sigmacoh*cohModel(x,params);
%     else
%         fcon = params.sigma*conModel(x,params);
%         fcoh = params.sigma*cohModel(x,params);
%     end
%     plot(x,fcon,'Color',clist(1,:));
%     plot(x,fcoh,'Color',clist(3,:));
%     % now plot the unattended curves    
%     title(sprintf('L: %2.3f.',likelihood));
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
    if pobs(12)==1 && pobs(8)==0 % left correct
        extra = params.left_correct;
    elseif pobs(12)==1 && pobs(8)==1 % right correct
        extra = params.right_correct;
    elseif pobs(12)==0 && pobs(8)==0
        extra = params.left_incorr;
    elseif pobs(12)==0 && pobs(8)==1
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

if isfield(params,'merge')
    cval = mean(abs(cons));
    mval = mean(abs(cohs));        
    pnoise = params.psigma*sqrt(abs(sum(abs(beta)*[cval ;mval])));
    
    if params.merge<0, params.merge=0; end
    if params.merge>1, params.merge=1; end
    
    tnoise = params.sigma * params.merge + pnoise * (1-params.merge);
    
    if obs(8)==1
        prob = normcdf(0,effect,tnoise,'upper');
    elseif obs(8)==0
        prob = normcdf(0,effect,tnoise);
    else warning('failure'); keyboard
    end
    
elseif params.poissonNoise
    
    cval = mean(abs(cons));
    mval = mean(abs(cohs));
        
    noise = sqrt(abs(sum(abs(beta)*[cval ;mval])));
    
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