
function r = cc_naka(i,params)
% Naka Rushton (based on 
% i: intensity
% sigma: semi-sat constant
% p: exponent
% q: exponent
% optional
% a: response gain (Def = 1)
% b: baseline shift (Def = 0)
% s: intensity gain (Def = 1)

if ~isfield(params,'sigma') || ~isfield(params,'N') || ~isfield(params,'p') || ~isfield(params,'q')
    help cc_naka
    return
end

N = params.N; sigma = params.sigma; p = params.p; q = params.q;

if ~isfield(params,'a')
    a = 1;
else
    a = params.a;
end

if ~isfield(params,'b')
    b = 0;
else
    b = params.b;
end

if ~isfield(params,'s')
    s = 1;
else
    s = params.s;
end

r  = a.*N.*(((i+b).^(p+q))./((i+b).^q+(sigma./s).^q));
