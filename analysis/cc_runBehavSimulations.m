function simData = cc_runBehavSimulations()

%% This is a simulation function to compare different possible models.

%% Model A

% Model A is the following structure:
% I(L,con) = Stimulus Contrast Intensity on Left Side
% S(L,con) = f(I) naka-rushton eqn
% D(R,con) = S(R,con)-S(L,con) diff function
% P(R) = f(a * D(con) + (1-a) * D(coh)) | a is attention to contrast, f is
% the weibull function

%% Naka-Ruhston

simData.nakap.con_c = 0:.01:1; % contrast
simData.nakap.coh_c = 0:.01:1; % coherence

simData.nakap.con_n = 1; % slope parameter of contrast
simData.nakap.coh_n = 1;

simData.nakap.con_c50 = .1; % inflexion point
simData.nakap.coh_c50 = .2;

simData.nakap.con_k = 1; % scale?
simData.nakap.coh_k = 1;

simData.nakaCon = nakarushton(simData.nakap.con_c,simData.nakap.con_n,simData.nakap.con_c50,simData.nakap.con_k);
simData.nakaCoh = nakarushton(simData.nakap.coh_c,simData.nakap.coh_n,simData.nakap.coh_c50,simData.nakap.coh_k);

figure
hold on
map = brewermap(15,'PuOr');
plot(simData.nakap.con_c,simData.nakaCon,'-','Color',map(5,:));
plot(simData.nakap.coh_c,simData.nakaCoh,'-','Color',map(11,:));

%% Generate Randomized Trials
trials = 1;

lCon = []; rCon = []; lCoh = []; lCon = []; task = [];
for t = 1:trials

%% Functions

% def: function w = weibull(x,params) params = T, Beta, lambda, min

function r = nakarushton(c,n,c50,k,cn)

% c is the contrast range
% n is the slope
% c50 is the inflexion point

if ieNotDefined('cn'),cn = c;end

r  = k*(c.^n./(cn.^n+c50.^n));