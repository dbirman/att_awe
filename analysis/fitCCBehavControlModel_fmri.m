function fit = fitCCBehavControlModel_fmri(adata,info,crossval)
% CCBehavModel
%
% Copy of fitCCBehavControlModel that allows you to fit a doublesigma (one
% sigma for contrast, one for coherence)

ROIs = {'V1','V2','V3','V4','V3a','V3b','V7','MT'};

global fixedParams

adata = adata(~any(isnan(adata),2),:);
osize = size(adata,1);

adata = adata(adata(:,9)==-1,:);
% disp(sprintf('Reducing data to %i control trials from %i',size(adata,1),osize));

%% Parse input

sigmaval = info.sigma;
model = info.model;
lapserate = info.lapse;

% if ROI model, pull out ROIs
if strfind(info.model,'roi')
    rois = ROIs(info.rois);
    
    roifit = struct;
    for ri = 1:length(rois)
        roifit.(rois{ri}) = struct;
        if strfind(info.model,'att')
            fixedParams.att = 1;
            for ai = 1:2
                roifit.(rois{ri}).confit(ai,:) = info.respcon(info.rois(ri),ai,:);
                roifit.(rois{ri}).cohfit(ai,:) = info.respcoh(info.rois(ri),ai,:);
            end
        else
            roifit.(rois{ri}).confit = info.respcon(info.rois(ri),:);
            roifit.(rois{ri}).cohfit = info.respcoh(info.rois(ri),:);
        end
    end
else
    confit = squeeze(mean(info.respcon(info.conGroup,:),1));
    cohfit = squeeze(mean(info.respcoh(info.cohGroup,:),1));
end

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
        
        trainfit = fitCCBehavControlModel_fmri(traindata,info,0);
%         trainr2(fold) = trainfit.pseudor2;
        sigmas(fold) = trainfit.params.sigma;
        testinfo = info;
        testinfo.fitmodel = trainfit;
        testfit = fitCCBehavControlModel_fmri(testdata,testinfo,0);
%         testr2(fold) = testfit.pseudor2;
        aresp = [aresp ;testdata(:,8)];
        aprobs = [aprobs testfit.probs];
        pseudor2(fold) = testfit.pseudor2;
        like(fold) = testfit.likelihood;
        nulllike(fold) = testfit.null.likelihood;
    end
    
    fit = fitCCBehavControlModel_fmri(adata,info,0);
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

%% Setup responses
fixedParams.x = 0:.001:1;
if strfind(model,'roi')
    % confit/cohfit contain the relevant data
    fixedParams.roifit = roifit;
    fixedParams.rois = rois;
    fixedParams.roi = length(rois);
else
    fixedParams.con = confit;
    fixedParams.coh = cohfit;
    fixedParams.roi = 0;
end
        
%% Special condition: just getting BIC for a model
if isfield(info,'fitmodel') && isstruct(info.fitmodel)
    fit = struct;
    
    nullinfo = info;
    nullinfo.fitmodel = [];
    if strfind(model,'roi')
        if strfind(model,'onebeta')
            nullinfo.model = 'null,roi,onebeta';
        else
            nullinfo.model = 'null,roi';
        end
    else
        nullinfo.model = 'null';
    end
    nullinfo.lapse = 0;
    if strfind(model,'onebeta')
        fit.null = struct;
        fit.null.likelihood = 0;
    else
        fit.null = fitCCBehavControlModel_fmri(adata,nullinfo,0);
    end
    [fit.likelihood,fitted] = fitBehavModel(info.fitmodel.params,adata,-1);
    probs = fitted.probs;
    fit.adata = fitted.adata;
    fit.BIC = 2*fit.likelihood + info.fitmodel.numParams * log(size(adata,1));
    fit.AIC = 2*fit.likelihood + info.fitmodel.numParams * 2;
    fit.params = info.fitmodel.params;
    fit.probs = probs;
    fit.resp = adata(:,8);
    fit.cd = nanmean(probs(fit.resp==1))-nanmean(1-probs(fit.resp==0));
    fit.pseudor2 = 1 - (fit.likelihood / fit.null.likelihood);
    return
end

%% Setup model parameters
if strfind(model,'null')
    % ONLY ALLOW BIAS TO CHANGE (intercept only model)
    initparams.conmodel = 4;
    initparams.cohmodel = 4;
    if strfind(model,'roi')
        if strfind(model,'onebeta')
            fixedParams.onebeta = 1;
            for ri = 1:length(rois)
                beta = sprintf('beta_control_%s_w',rois{ri});
                
                initparams.(beta) = 0;
            end
        else
            for ri = 1:length(rois)
                cbeta = sprintf('beta_control_%s_conw',rois{ri});
                mbeta = sprintf('beta_control_%s_cohw',rois{ri});
                initparams.(cbeta) = 0;
                initparams.(mbeta) = 0;
            end
        end
    else
        initparams.beta_control_con_conw = 0;
        initparams.beta_control_con_cohw = 0;%[0 -1 1];
        initparams.beta_control_coh_cohw = 0;
        initparams.beta_control_coh_conw = 0;%[0 -1 1];
    end
    initparams.bias = [0 -1 1];
    initparams.sigma = 1;
    initparams.coh_gain = 0;
    initparams.poissonNoise = 0;
    initparams.lapse = 0;
    [~, fit] = fitModel(initparams,adata,-1);
    return
elseif strfind(model,'freeze')
    % ONLY ALLOW SIGMA TO CHANGE

    initparams.conmodel = 4;
    initparams.cohmodel = 4;
    if strfind(model,'roi')
        warning('NOT IMPLEMENTED');
        for ri = 1:length(rois)
            cbeta = sprintf('beta_control_%s_conw',rois{ri});
            mbeta = sprintf('beta_control_%s_cohw',rois{ri});
            initparams.(cbeta) = rand;
            initparams.(mbeta) = rand;
        end
    else
        initparams.beta_control_con_conw = 1;
        initparams.beta_control_con_cohw = 0;%[0 -1 1];
        initparams.beta_control_coh_cohw = 1;
        initparams.beta_control_coh_conw = 0;%[0 -1 1];
    end
    initparams.bias = 0;%;[0 -1 1];
    initparams.sigma = [0.1 eps 1];
    initparams.poissonNoise = 0;
    initparams.lapse = lapserate;
    [~, fit] = fitModel(initparams,adata,-1);
    return
elseif strfind(model,'sigma')
%     disp('Fitting sigma...');
    % SPECIAL CONDITION: Fitting sigma parameter
    initparams.conmodel = 4;
    initparams.cohmodel = 4;
    if strfind(model,'roi')
        if strfind(model,'onebeta')
            fixedParams.onebeta = 1;
            for ri = 1:length(rois)
                beta = sprintf('beta_control_%s_w',rois{ri});
                
                initparams.(beta) = [rand -inf inf];
            end
        else
            for ri = 1:length(rois)
                cbeta = sprintf('beta_control_%s_conw',rois{ri});
                mbeta = sprintf('beta_control_%s_cohw',rois{ri});
    %             if strfind(rois{ri},'V1')
    %                 initparams.(cbeta) = 1;
    %             else
                initparams.(cbeta) = [rand -inf inf];
    %             end
                initparams.(mbeta) = [rand -inf inf];
            end
        end
    else
        initparams.beta_control_con_conw = 1;
        initparams.beta_control_con_cohw = [0 -1 1];
        initparams.beta_control_coh_cohw = 1;
        initparams.beta_control_coh_conw = [0 -1 1];
    end
    initparams.bias = [0 -1 1];
    if strfind(model,'gain')
        initparams.coh_gain = [1.2 0 inf];
        initparams.sigma = sigmaval;
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
                if strfind(model,'roi')
%                     warning('SIGMA PARAMETER NOT INTERPRETABLE: FROZEN AT 1');
                    initparams.sigma = 1;
                else
                    initparams.sigma = [sigmaval eps 1];
                end
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
                if strfind(model,'roi')
%                     warning('SIGMA PARAMETER NOT INTERPRETABLE: FROZEN AT 1');
                    initparams.sigma = 1;
                else
                    initparams.sigma = [sigmaval eps 1];
                end
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
[fit.likelihood, fitted] = fitBehavModel(bestparams,adata,0);
probs = fitted.probs;
fit.adata = fitted.adata;
fit.BIC = 2*fit.likelihood + length(bestparams) * log(size(adata,1));
fit.AIC = 2*fit.likelihood + length(bestparams) * 2;
fit.probs = probs;
fit.resp = adata(:,8);
fit.cd = nanmean(probs(fit.resp==1))-nanmean(1-probs(fit.resp==0));
fit.numParams = length(bestparams);

function [likelihood, fit] = fitBehavModel(params,adata,f)
%%
global fixedParams

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

if fixedParams.roi
    betas = zeros(2,fixedParams.roi);
    if fixedParams.onebeta
        for ri = 1:fixedParams.roi
            betas(:,ri) = repmat(params.(sprintf('beta_control_%s_w',fixedParams.rois{ri})),2,1);
        end
    else
        for ri = 1:fixedParams.roi
            betas(1,ri) = params.(sprintf('beta_control_%s_cohw',fixedParams.rois{ri}));
            betas(2,ri) = params.(sprintf('beta_control_%s_conw',fixedParams.rois{ri}));
        end
    end
else
    % compute betas
    betas = zeros(6,2);
    betas(1,:) = [params.beta_control_coh_conw params.beta_control_coh_cohw];
    betas(2,:) = [params.beta_control_con_conw params.beta_control_con_cohw];
    % betas(3,:) = [params.beta_att_coh_conw params.beta_att_coh_cohw];
    % betas(4,:) = [params.beta_unatt_coh_conw params.beta_unatt_coh_cohw];
    % betas(5,:) = [params.beta_att_con_conw params.beta_att_con_cohw];
    % betas(6,:) = [params.beta_unatt_con_conw params.beta_unatt_con_cohw];
end

% compute effects
if fixedParams.roi
    % compute the effect for each ROI for each side
    roiEff = zeros(size(adata,1),fixedParams.roi);
    for ri = 1:fixedParams.roi
        if fixedParams.att
            % ATTENTION MODEL 
            
            % index response by whether attention is directed to contrast
            % or coherence
            for i = 1:size(adata,1)
                fixedParams.con = fixedParams.roifit.(fixedParams.rois{ri}).confit(adata(i,1),:);
                fixedParams.coh = fixedParams.roifit.(fixedParams.rois{ri}).cohfit(adata(i,1),:);
                
                conEffL = (conModel(adata(i,4),params)-conModel(adata(i,2),params));
                conEffR = (conModel(adata(i,5),params)-conModel(adata(i,2),params));
                conEff = conEffR - conEffL;
                cohEffL = (cohModel(adata(i,6),params)-cohModel(adata(i,3),params));
                cohEffR = (cohModel(adata(i,7),params)-cohModel(adata(i,3),params));
                cohEff = cohEffR - cohEffL;

                roiEff(i,ri) = conEffR + cohEffR - conEffL -cohEffL;
            end
        else
            fixedParams.con = fixedParams.roifit.(fixedParams.rois{ri}).confit;
            fixedParams.coh = fixedParams.roifit.(fixedParams.rois{ri}).cohfit;
            conEffL = (conModel(adata(:,4),params)-conModel(adata(:,2),params));
            conEffR = (conModel(adata(:,5),params)-conModel(adata(:,2),params));
            conEff = conEffR - conEffL;
            cohEffL = (cohModel(adata(:,6),params)-cohModel(adata(:,3),params));
            cohEffR = (cohModel(adata(:,7),params)-cohModel(adata(:,3),params));
            cohEff = cohEffR - cohEffL;

            roiEff(:,ri) = conEffR + cohEffR - conEffL -cohEffL;
        end
    end
    
else
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
    roiEff = zeros(size(conEff));
end

probs = zeros(1,size(adata,1));

for ai = 1:size(adata,1)
    obs = adata(ai,:);
    
    if length(conEff)>1
        ce = conEff(ai);
        me = cohEff(ai);
        lr_con = [conEffL(ai) conEffR(ai)];
        lr_coh = [cohEffL(ai) cohEffR(ai)];
    else
        ce = []; me = [];
        lr_con = []; lr_coh = [];
    end
    
    if ai>1
        prob = getObsProb(obs,params,adata(ai-1,:),betas,ce,me,lr_con,lr_coh,roiEff(ai,:));
    else
        prob = getObsProb(obs,params,[],betas,ce,me,lr_con,lr_coh,roiEff(ai,:));
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

fit.probs = probs;
fit.adata = adata;

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

function prob = getObsProb(obs,params,pobs,betas,conEff,cohEff,cons,cohs,roiEff)
%%
global fixedParams
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

if fixedParams.roi
    effect = beta * roiEff' + params.bias + extra;
else
    effect = beta * [conEff cohEff]' + params.bias + extra;
end

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
    if fixedParams.roi
        noise = sqrt(abs(sum(abs(beta)*roiEff')));
    else
        cval = mean(abs(cons));
        mval = mean(abs(cohs));

        noise = sqrt(abs(sum(abs(beta)*[cval ;mval])));
    end
    
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