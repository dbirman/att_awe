function matlab_full()
%MATLAB_FULL Run the Full Simulation
%   This is a simulation function to compare different possible models. It
%   will load the neural and behavioral datasets stored in fullData.mat

%% Load and Initialize Parameters

% Naka-Rushtonish: N * i^(p+q) / [i^q * (sigma/s)^q]
params = struct; % defaults
% nake def
params.def.naka.p = 0.1; params.def.naka.q = 0.1; params.def.naka.N = 1.5; params.def.naka.sigma = 0.5;
params.def.linear.b0 = 1; params.def.linear.a = 0.5; 
% params.con = params.def.naka;
% params.coh = params.def.naka;
params.con = params.def.linear;
params.coh = params.def.linear;

global settings
settings.flag = 'linear';

% lb = [0,0,0,0,1,0,1]; lb = [lb lb];
% ub = [10,10,10,10,1.01,0.01,1.01]; ub = [ub ub];
lb = [0,-5];lb = [lb lb];
ub = [5,5]; ub = [ub ub];

fd = loadFullData;

tasks = {'coherence','contrast'};

fd.rois = {'V3a','V3b','V1','V2','V3','MT'};

%% Step 1: Compute Parameter Fits for Linear Eqn by Subj, ROI, and Task
params = struct;
roiparams = struct;
tasks = {'coherence','contrast'};
% group by subject
for si = 1:length(fd.nsubjects)
    % group by ROI
    for ri = 1:length(fd.rois)
        % group by task
        for ti = 1:2
            % select
            data = select(fd.ndata,7,mrStr2num(strrep(fd.nsubjects{si},'s','')));
            data = select(data,6,ri);
            data = select(data,3,ti);           
            
            % data is amp, con, coh
            comb = data(:,[1 4 5]);
            
            X = [ones(size(comb,1),1) comb(:,2:3)]; Y = comb(:,1);
            
            bL = X\Y;
            
            params.con.b0 = bL(1); params.coh.b0 = bL(1);
            params.con.a = bL(2); params.coh.a = bL(3);      
            
            roiparams.(strcat('s',num2str(mrStr2num(strrep(fd.nsubjects{si},'s',''))))).(tasks{ti}).(fd.rois{ri}) = params;
        end
    end
end

%% Now Pull the Behavioral Data out by group (we need resp, con, coh)
adata = struct;
for si = 1:length(fd.subjects)
    sid = fd.subjects{si};
    data = select(fd.data,14,mrStr2num(strrep(sid,'s','')));
    data = select(data,8,1);
    
    % now we pull out the data for each subgroup
    % dat.(task).(condition)
    ccon = data(:,5)-data(:,4);
    ccoh = data(:,7)-data(:,6);
    ccatch = data(:,11);
    resp = data(:,12);
    nocatch = data(:,13);
    side = data(:,15);
    task = data(:,3);
    
    nocatch = logical(nocatch); catchs = ccatch==1;
    
    runs = {'nocatch','catch'};
    tasks = {'coherence','contrast'};
    catchs = {'cued','miscued'};
    
    dat = struct;
    for ri = 1:2
        for ti = 1:2
            if ri==1
                selector = logical((ccatch==-1).*(task==ti));
                dat.(tasks{ti}).nocatch.con = ccon(selector);
                dat.(tasks{ti}).nocatch.coh = ccoh(selector);
                dat.(tasks{ti}).nocatch.resp = resp(selector);
                dat.(tasks{ti}).nocatch.side = side(selector);
            else
                for ci = 0:1
                    selector = logical((task==ti).*(ccatch==ci));
                    dat.(tasks{ti}).(catchs{ci+1}).con = ccon(selector);
                    dat.(tasks{ti}).(catchs{ci+1}).coh = ccoh(selector);
                    dat.(tasks{ti}).(catchs{ci+1}).resp = resp(selector);
                    dat.(tasks{ti}).(catchs{ci+1}).side = side(selector);
                end
            end
        end
    end
    
    % data has no been pulled
    adata.(sid) = dat;
end

%% Move the ROI Data into bdata

for si = 1:length(fd.subjects) % only for subjs in behav
    sid = fd.subjects{si};
    
    for ti = 1:2
        
        task = tasks{ti};
        
        copyFields = fields(roiparams.(sid).(task));
        
        for cfi = 1:length(copyFields)
            adata.(sid).(task).(copyFields{cfi}) = roiparams.(sid).(task).(copyFields{cfi});
        end
    end
end

%% Go through the data and get the relevant V1/MT intercept + slope

for si = 1:length(fd.subjects) % only for subjs in behav
    sid = fd.subjects{si};
    
    for ti = 1:length(tasks)
        task = tasks{ti};
        
        % run the task effect estimator
        adata.(sid).(task) = cc_fullEstimator(adata.(sid).(task));
    end
end

stop = 1; 

function data = cc_fullEstimator(data)
%% setup
% estimate the effect size in each condition, e.g. the slope of the d' line
conds = {'nocatch','cued','miscued'};
group = [1 1 2];
groups = {[1 2], [2 1]};
rois = {'V3a','V3b','V1','V2','V3','MT'};

%% run each cond
for ci = 1:length(conds)
    cond = conds{ci};
    
    dat = data.(cond);
    
    pars = cc_fitCumGauss(dat.con',dat.coh',dat.resp',0);
    
%     con_range = -0.3:.01:.3;
%     con_eff = norminv(normcdf(con_range,pars(1),pars(2)));
%     con_range = con_range(logical((con_eff~=Inf).*(con_eff~=-Inf)));
%     con_eff = con_eff(logical((con_eff~=Inf).*(con_eff~=-Inf)));

    dat.conEff = 1/pars(2);
    dat.cohEff = 1/pars(4);
    
    %% responses
    conResp = [];
    cohResp = [];
    conName = {};
    cohName = {};
    
    for ri = 1:length(rois)
        roi = rois{ri};
        
        % get all the responses
        conResp = [conResp data.(roi).con.a];
        conName{end+1} = strcat(roi,'_con');
        cohResp = [cohResp data.(roi).coh.a];
        cohName{end+1} = strcat(roi,'_coh');
    end
    
    %% estimate noise in each region if it were solely responsible
    conNoise = conResp / dat.conEff;
    cohNoise = cohResp / dat.cohEff;
    
    %% estimate noise in each region if decoding was linear
    conNoiseShared = conResp ./ (dat.conEff * conResp ./ sum(conResp));
    cohNoiseShared = cohResp ./ (dat.cohEff * cohResp ./ sum(cohResp));
    
    dat.conNoise = conNoise;
    dat.cohNoise = cohNoise;
    dat.conNoiseShared = conNoiseShared;
    dat.cohNoiseShared = cohNoiseShared;
    %% return data
    data.(cond) = dat;
end


function f2 = plotResp2D(con_range,coh_range,params)

% we are going to plot two plots, each will be the average response 
[x,y,out] = respRange(con_range,coh_range,params);

% columns are contrast, rows are coherence

f2 = figure;
subplot(2,1,1) % contrast
plot(x,mean(out,1));
title('contrast')
subplot(2,1,2)
plot(y,mean(out,2));
title('coherence')


function resid = respResid(dat,v)
Rt = dat(:,1); con = dat(:,2); coh = dat(:,3);

resid = (respModel(con,coh,v) - Rt).^2;

function R = respModel(con,coh,v)
global settings

if length(con) > 1
    R = arrayfun(@(x1,x2) respModel(x1,x2,v), con, coh);
else
    if strcmp(settings.flag,'naka')
        R = respNaka(con,v(1),v(2),v(3),v(4),v(5),v(6),v(7)) + ...
            respNaka(coh,v(8),v(9),v(10),v(11),v(12),v(13),v(14));
    elseif strcmp(settings.flag,'linear')
        R = respLinear(con,coh,v(1),v(2),v(4));
    end
end

function f3 = plotResp3D(con_range,coh_range,params)

[~,~,out] = respRange(con_range,coh_range,params);
f3 = figure;
surf(coh_range,con_range,out)
colorbar
ylabel('Coherence')
xlabel('Contrast')

function [x,y,out] = respRange(con_range,coh_range,params)
coh_range = fliplr(coh_range);
out = zeros(length(con_range),length(coh_range));
v = params2vec(params);
for i=1:length(coh_range)
    for j=1:length(con_range)
        out(i,j) = respModel(con_range(j),coh_range(i),v);
    end
end
x = con_range;
y = coh_range;

function y = respLinear(iCon,iCoh,b0,aCon,aCoh)
y = b0 + aCon*iCon + aCoh*iCoh;

function y = respNaka(i,p,q,N,sigma,a,b,s)
% p/q exponents
% N scale
% sigma semi-saturation constant
% a response gain def=1
% b baseline shift def=0
% s intensity gain def=1
y = a * N * (((i+b)^(p+q)) / (((i+b)^q)+((sigma/s)^q)));

function data = select(data,col,value)

scol = data(:,col);
data = data(scol==value,:);

function params = vec2params(v)
global settings
if strcmp(settings.flag,'naka')
    params.con.p = v(1);
    params.con.q = v(2);
    params.con.N = v(3);
    params.con.sigma = v(4);
    params.con.a = v(5);
    params.con.b = v(6);
    params.con.s = v(7);
    params.coh.p = v(8);
    params.coh.q = v(9);
    params.coh.N = v(10);
    params.coh.sigma = v(11);
    params.coh.a = v(12);
    params.coh.b = v(13);
    params.coh.s = v(14);
elseif strcmp(settings.flag,'linear')
    params.con.b0 = v(1);
    params.con.a = v(2);
    params.coh.b0 = v(3);
    params.coh.a = v(4);
end

function v = params2vec(params)
global settings

if strcmp(settings.flag,'naka')
    if isfield(params.con,'a'), a1 = params.con.a; else a1 = 1; end
    if isfield(params.con,'b'), b1 = params.con.b; else b1 = 0; end
    if isfield(params.con,'s'), s1 = params.con.s; else s1 = 1; end
    if isfield(params.coh,'a'), a2 = params.coh.a; else a2 = 1; end
    if isfield(params.coh,'b'), b2 = params.coh.b; else b2 = 0; end
    if isfield(params.coh,'s'), s2 = params.coh.s; else s2 = 1; end

    v = [params.con.p,params.con.q,params.con.N,params.con.sigma,a1,b1,s1,...
        params.coh.p,params.coh.q,params.coh.N,params.coh.sigma,a2,b2,s2];
elseif strcmp(settings.flag,'linear')
    v = [params.con.b0,params.con.a,params.coh.b0,params.coh.a];
else
    disp('unsure');
end
