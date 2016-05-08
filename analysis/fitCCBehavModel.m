function [params, fit] = fitCCBehavModel(adata,figs)
% CCBehavModel
%
% Fit the contrast (naka-rushton) and coherence (linear) models to the data
% obtained from the behavioral experiment. Roughly we will do the
% following:
%
% Figure out the effect size:
% 
%   Rmax
%   c50
%   n
%
%   slope
%
%   sigma_con <- noise in the contrast response
%   sigma_coh <- noise in the coherence response

% remove any NaN
adata = adata(~any(isnan(adata),2),:);

%% Parameters
global fixedParams
fixedParams = struct;

initparams.Rmax = [1 -inf inf];
initparams.c50 = [0.5 0 1];
initparams.n = 1;

initparams.slope = [1 -inf inf];

initparams.sigma = 1;

initparams.alphacon = [0.9 0 1];
initparams.alphacoh = [0.9 0 1];
initparams.alphacon_un = [0.8 0 1];
initparams.alphacoh_un = [0.8 0 1];

% initparams.unscalecon = [1 0 1];
% initparams.unscalecoh = [1 0 1];

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

[params, fit] = fitModel(initparams,adata,f);

function [bestparams,fit] = fitModel(params,adata,f)

[initparams, minparams, maxparams] = initParams(params);

bestparams = fmincon(@(p) fitBehavModel(p,adata,f),initparams,[],[],[],[],minparams,maxparams);

fit.params = getParams(bestparams);
[fit.likelihood, fit.fit] = fitBehavModel(params,adata,0);

function likelihood = fitBehavModel(params,adata,f)
%%
params = getParams(params);

% validate params
if params.sigma < 0, params.sigma = 0; end
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
    title(sprintf('c50: %0.2f slope: %1.2f alphacon %0.2f alphacoh %0.2f',params.c50,params.slope,params.alphacon,params.alphacoh));
end

for ai = 1:size(adata,1)
    
    prob = -1;
    obs = adata(ai,:);
    
    conEff = (conModel(obs(5),params)-conModel(obs(2),params)) - (conModel(obs(4),params)-conModel(obs(2),params));
    cohEff = (cohModel(obs(7),params)-cohModel(obs(3),params)) - (cohModel(obs(6),params)-cohModel(obs(3),params));
    conProb = normcdf(0,conEff,params.sigma);
    cohProb = normcdf(0,cohEff,params.sigma);
    if obs(8)
        conProb = 1-conProb;
        cohProb = 1-cohProb;
    end
    
    switch obs(1)
        case 1
            %COHERENCE CONTROL
            prob = cohProb * params.alphacoh + conProb * (1-params.alphacoh);
        case 2
            %CONTRAST CONTROL
            prob = conProb * params.alphacon + cohProb * (1-params.alphacon);
        case -1
            % COHERENCE MAIN
            if obs(9)==0
                % main
                prob = cohProb * params.alphacoh_un + conProb * (1-params.alphacoh_un);
            elseif obs(9)==1
                % catch
                
            end
        case -2
            % CONTRAST MAIN
            if obs(9)==0
                % main
                prob = conProb * params.alphacon_un + cohProb * (1-params.alphacon_un);
            elseif obs(9)==1
                % catch
            end
    end
    
    probs(ai) = prob;
    if prob >= 0
        likelihood = likelihood + log(prob);
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
end

function out = conModel(con,params)

out = params.Rmax .* ((con.^params.n) ./ (con.^params.n + params.c50.^params.n)); 

function out = cohModel(coh,params)

out = params.slope .* coh;

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