function findbestgroupmodel( control,attend )

% we'll use maximum likelihood to find the best fitting contrast and
% coherence models, we'll compare various exponents and linear models, but
% keep in mind this is all fit to the threshold data, so it's an
% approximation at best.

% not clear whether to do this on the average or on all the data? Probably
% better on all the data...

%% Fit Contrast

basecon = 0.25;
con = [0.325 0.4 0.55 0.85 0.4];
stim.con = con;
stim.basecon = basecon;
data = [squeeze(control(:,2,:)),attend(:,2)];
confits = tryModels(stim,data);

%% Fit Coherence Functoin
basecoh = 0;
coh = [0.15 0.3 0.45 0.6 0.3];
data = [squeeze(control(:,1,:)),attend(:,1)];
stim.con = coh; stim.basecon = basecoh;
cohfits = tryModels(stim,data);

%%
h = figure;
subplot(211), hold on % contrast
cmap = brewermap(3,'PuOr');
plot(0:.01:1,repmat(confits.int,1,101),'--','Color',cmap(1,:));
plot(confits.con,confits.data,'o','MarkerFaceColor',cmap(1,:),'MarkerEdgeColor',[1 1 1]);
plot(0:.01:1,(0:.01:1)*confits.b(2)+confits.b(1),'-','Color',cmap(1,:));
axis([0 1 0 0.2])
xlabel('Base Contrast')
ylabel('Just Noticeable Difference')
[r,p] = corr(confits.con,confits.data);
title(sprintf('Contrast: R^2 %2.1f%% p = %0.3f',r^2*100,p));
drawPublishAxis
subplot(212), hold on
plot(cohfits.con,cohfits.data,'o','MarkerFaceColor',cmap(3,:),'MarkerEdgeColor',[1 1 1]);
plot(0:.01:1,repmat(cohfits.int,1,101),'--','Color',cmap(3,:));
plot(0:.01:1,(0:.01:1)*cohfits.b(2)+cohfits.b(1),'-','Color',cmap(3,:));
axis([0 1 0 0.4])
xlabel('Base Coherence')
ylabel('Just Noticeable Difference')
[r,p] = corr(cohfits.con,cohfits.data);
title(sprintf('Contrast: R^2 %2.1f%% p = %0.3f',r^2*100,p));
drawPublishAxis

fname = fullfile(datafolder,'avg_groupmodel.pdf');
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'InvertHardCopy','off');
set(gcf,'Color',[1 1 1]);
set(gca,'Color',[1 1 1]);
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fname,'-dpdf');

function fits = tryModels(lstim,data)
global stim
stim = lstim;
con = repmat(lstim.con,size(data,1),1);

con = con(:);
data = data(:);
con = con(~isnan(data));
data = data(~isnan(data));
% fit the intercept model
intercept = nanmean(data(:));
slope = 1/intercept; % this is the slope of the response model needed for those thresholds
ss1 = sum((data(:)-intercept).^2);
c2 = sum((data(:)-intercept).^2./intercept);
% fit the linear modeal
data = data(:); 
con = con(:);
con = [ones(size(con)) con];
b = con\data; % b(1) = intercept, b(2) = slope
slope2 = (0:.01:1)*b(2)+b(1);
estimate = con(:,2)*b(2)+b(1);
ss2 = sum((data-estimate).^2);
c22 = sum((data-estimate).^2./estimate);

p = 1-chi2pdf(c2-c22,1);

fits.p = p;
fits.int = intercept;
fits.b = b;
fits.data = data;
fits.con = con(:,2);
fits.ss = [ss1 ss2];
% initparams.conslope = [1 -inf inf];
% initparams.conmodel = 1;
% params = initParams(initparams);
% 
% optimParams = optimset('Algorithm','levenberg-marquardt','MaxIter',inf,'Display','off','DiffMinChange',0.1);
% 
% linparams = lsqnonlin(@(p) conResidual(p,data),params,-inf(size(params)),inf(size(params)),optimParams);
% linres = conResidual(linparams,data);
% linearfit.params = getParams(linparams);
% linearfit.res = linres;
% linearfit.n = 1; % num of params

%%
% nakafit = cell(1,3);
% for ni = 1:3
%     clear initparams
%     initparams.conRmax = [30 -inf inf];
%     initparams.conc50 = [0.75 0 1];
%     initparams.conn = ni;
%     initparams.conmodel = 2;
%     params = initParams(initparams);
%     nakaparams = lsqnonlin(@(x) conResidual(x,data),params,-inf(size(params)),inf(size(params)),optimParams);
%     nakares = conResidual(nakaparams,data);
%     nakafit{ni}.params = getParams(nakaparams);
%     nakafit{ni}.res = nakares;
%     nakafit{ni}.n = 3; % num of params
% end
% 
% % fits.lin = linearfit;
% fits.naka = nakafit;

function res = conResidual(params,data)
%%
global stim

params = getParams(params);

x = -5:.01:5;
r = conModel(x,params);

estimate = zeros(size(stim.con));
for ci = 1:length(stim.con)
    con = stim.con(ci);
    idx = find(x>=con,1);
    rcon = r(idx);
    rthresh = rcon+1; % sigma is 1, so d'=1 is +1 on the resp
    tidx = find(r>=rthresh,1);
    estimate(ci) = x(tidx);
end

res = data-repmat(estimate,size(data,1),1);
res = res(:);
res(isnan(data(:)))=0;

if any(isnan(res)), keyboard; end
% function likelihood = cohResidual(params,data)

function [initparams, minparams, maxparams] = initParams(params)

global fixedParams
fixedParams = struct;

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