% fitgamma.m
%
%      usage: fitgamma(hdr{},<conv_timing>,<time>,<minfit>,<maxfit>,<dispfit>)
%         by: justin gardner, n by dan birman
%       date: 08/21/03
%    purpose: fit a gamma function to hemodynamic response
%             the parameters returned are
%             [amp tau timelag offset exponent]
%             offset is usually forced to zero
%             exponent is either 5 or 6
%             minfit is an array of length three specifying
%             the minimum values for amp tau and timelag respectively
%             maxfit is the maximum
%
%        e.g.:
%tr = 0.8;
%time = tr:tr:25*tr;
%hrf = mygamma(time-3,5,1);hrf = hrf/max(hrf);
%hrfnoise = hrf + rand(1,length(hrf))/2;
%bestfit = fitgamma(hrfnoise,time,[],[],1);
%
% This is a generalization of the fitgamma function that allows you to fit
% a single gamma function to many outputs. You input a set of HDR functions
% and a set of convolution filters (stimulus timing) that the HDR functions
% correspond to in theory. If you just want to fit multiple gamma functions
% at once just leave conv_timing empty.
%
function bestfit = fitgamman(hdr,conv_timing,time,minfit,maxfit,dispfit)

% check arguments
if ~any(nargin == [1 2 3 4 5])
  help fitgamma
  return
end

if ~exist('conv_timing','var') || isempty(conv_timing)
    conv_timing = [];
end
if ~exist('time','var') || isempty(time)
  time = 0:(length(hdr)-1);
end
if ~exist('minfit','var')
  minfit = [0.01   -10    0    0   -1  0];
end
if ~exist('maxfit','var')
  maxfit = [1  10    5       inf  1  50];
end
if ~exist('dispfit','var')
  dispfit = 0;
end
% whether to display errors or not
if (dispfit),displsqnonlin = 'final'; else displsqnonlin = 'off'; end

% make sure we have a column vector
if (size(hdr,1) == 1)
  hdr = hdr';
end

global numcalls;numcalls = 0;
maxiter = inf;
bestresnorm = inf;

% these will be set in the fits (not allowed to fix)
% so there value is passed to the objective function
% as these globals
global exponent;
global offset;

for i = 4:7
  % set fixed parameters
  exponent = i;
  offset = 0;
  % set initparams
  %             amp         tau lag alpha   lambda  beta
  initparams = [0.1         0.5   1   1.1       0       10];
  % fit function using lsqnonlin in LevenbergMarquardt mode.
  [fitparams resnorm residual exitflag output lambda jacobian] = lsqnonlin(@gammaerr,initparams,minfit,maxfit,optimset('Algorithm','levenberg-marquardt','MaxIter',maxiter,'Display',displsqnonlin),hdr,conv_timing,time,dispfit);

  % Taken from Numerical Recipies, 
  % the leastsq function seems to return the transposed gradient
  % instead of the jacobian...
  jacobian = jacobian'*jacobian;
  reducedChiSquared = (residual*residual')/(length(hdr)-length(initparams));
  covar = sqrt(reducedChiSquared * inv(jacobian));

  % if we have the best fit then keep it.
  if (resnorm < bestresnorm)
    bestfit.params = [fitparams offset exponent];
    bestfit.covar = covar;
    bestfit.output = output;
    bestresnorm = resnorm;
  end
end

[bestfit.err, bestfit.fit] = gammaerr(bestfit.params,hdr,conv_timing,time,dispfit);
bestfit.err = reshape(bestfit.err,size(hdr));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function for fittnig gamma to estimated hdrs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [err, fit] = gammaerr(fitparams,ehdr,conv_timing,time,dispfit)

% number of times this routine has been called
global numcalls;
numcalls = numcalls + 1;

% if we are passed three parameters then the
% rest come from globals
global exponent;
global offset;

amp = fitparams(1);
tau = fitparams(2);
timelag = fitparams(3);
alpha = fitparams(4);
lambda = fitparams(5);
beta = fitparams(6);
% else all parameters passed in

% calculate function
gammafun = mygamma(time-timelag,exponent,tau);
gammafun = (gammafun-min(gammafun)) ./ (max(gammafun)-min(gammafun));
fitfun = (amp*gammafun+offset);

% use the base function to compute all of the computed functions
err = zeros(size(ehdr));

usealpha = [0 0 0 0 0 1 1 1 1 1];
fitout = [];
for ei = 1:size(ehdr,1)
    
    % compute the base function, which is some scaling of the impulse
    % response
    basefun = fitfun * beta;
    
    % fitfun should be scaled in the 100% coherence condition
    if usealpha(ei)
        fitfun = fitfun * alpha;
    end
    
    % compute compfun which will be the output for this particular timing
    % set
    if isempty(conv_timing)
        % just the basefun + fitfun
        compfun = basefun + fitfun;
    else
        % build up the exponential decay convolution
        timing = conv_timing(ei,:);
        ftiming = timing(timing>0);
        ftiming = exp(-(ftiming)*lambda);
        timing(timing>0) = ftiming;
        % convolve to get the additive response due to the stimulus
        compfun_add = conv(fitfun,timing);
        % add it up
        compfun = basefun + compfun_add(1:size(fitfun,2));
    end
    
    err(ei,:) = ehdr(ei,:) - compfun;
    fitout = [fitout ; compfun];
end
err = err(:)';

fit.canonical = fitfun;
fit.base = fitfun*beta;
fit.add = compfun_add;
fit.fits = fitout;

ssres = sum(err.^2);
sstot = sum((ehdr(:)-mean(ehdr(:))).^2);
fit.r2 = 1- ssres/sstot;

if (dispfit && (10*(floor((numcalls-1)/10))==(numcalls-1)))
clf;
  plot(time,ehdr(:),'k.');
  hold on
  plot(time,fitfun,'r-');

  hline(0);
  yaxis(min(-1,-2.5*mean(std(ehdr'))),max(1.5,2.5*mean(std(ehdr'))));
  title(sprintf('amp=%0.2f tau=%0.2f lag=%0.2f funcalls=%i',fitparams(1),fitparams(2),fitparams(3),numcalls));
  xlabel('time (sec)');
  ylabel('BOLD response (%)');
  drawnow
end