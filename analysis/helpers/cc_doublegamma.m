
%%
function out = cc_doublegamma(time,params)

% global fixedParams
% 
% if fixedParams.diff
    n = params.exponent;
    tau1 = params.tau1;
    amp1 = params.amp1;
    tau2 = params.tau2;
    amp2 = params.amp2;
    
    time1 = time-params.timelag1;
    out1 = ((time1/tau1).^(n-1).*exp(-time1/tau1))./(tau2*factorial(n-1));
    out1(time1<0) = 0;
    out1 = (out1-min(out1))./ (max(out1)-min(out1));
    out1 = amp1*out1;
    
    time2 = time-params.timelag2;
    out2 = ((time2/tau2).^(n-1).*exp(-time2/tau2))./(tau2*factorial(n-1));
    out2(time2<0) = 0;
    out2 = (out2-min(out2))./(max(out2)-min(out2));
    out2 = amp2*out2;
    
    out = out1+out2;
% else
%     n = params.exponent;
%     tau = params.tau;
%     time = time-params.timelag;
%     out = ((time/tau).^(n-1).*exp(-time/tau))./(tau*factorial(n-1));
%     out(time < 0) = 0;
%     out = (out-min(out)) ./ (max(out)-min(out));
%     out = params.amplitude*out;
% end

out = out/sum(out)/2; % normalize to sum=1 for 1% signal change / s
