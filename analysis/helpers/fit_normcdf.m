function fit = fit_normcdf( x, y )

guess = [0 1 0];

out = lsqnonlin(@(p) res_normcdf(p,x,y),guess);

figure(1);
mu = out(1);
sigma = out(2);
lapse = out(3);

[rres,fit] = res_normcdf(out,x,y);

params.mu = mu; params.sigma = sigma; params.lapse = lapse;
fit.params = params;
fit.sse_tot = sum((y-0.5).^2);
fit.r2 = 1 - sum(rres.^2)/fit.sse_tot;

function [res,fit] = res_normcdf(params, x,y)

fit.x = x;
fit.y = y;

if params(2)<=0, res = Inf*ones(size(y)); return; end

y_ = (1-2*params(3)) * normcdf(x,params(1),params(2)) + params(3);

fit.y_ = y_;

res = y_-y;

figure(1); clf; hold on

title(sum(res));
plot(x,y,'ok');
plot(x,y_,'--k');
