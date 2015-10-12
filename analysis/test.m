function test(con,coh,resp)
%% testing

nans = (isnan(resp)+isnan(con)+isnan(coh))>0;
con = con(~nans);
coh = coh(~nans);
resp = resp(~nans);

%% 
            cone = .3; cohe = .8; cont = 0.025; coht = 0.1;
            conbins = -cone:cont:cone; conrange = -cone+cont/2:cont:cone-cont/2;
            cohbins = -cohe:coht:cohe; cohrange = -cohe+coht/2:coht:cohe-coht/2;

            %%
r = corr(con',coh');
t = r / sqrt((1-r^2)/(length(con)-2));
p = 1-tcdf(t,length(con)-2); 
disp(sprintf('Contrast and Coherence r = %0.2f, p = %0.2f',corr(con',coh'),p));

if p > 0.05
    [cons, i] = sort(con);
    resps = resp(i);
    respc = cumsum(resps)/sum(resps);
    conBinned = binData(resp,con,conbins);
    conMean = cellfun(@nanmean,conBinned);
    conN = cellfun(@length,conBinned);
    cons = conbins(conN>1);
    respc = conMean(conN>1);
    conN = conN(conN>1);
    figure, hold on
    plot(cons,respc,'*b')
    axis([-.3 .3 0 1]);
    vline(0)
    hline(0.5)
    y = normcdf(-0.3:0.01:.3,mean(cons),std(cons));
    plot(-0.3:.01:.3,y,'r');
    
    % Initial Guesses
%     a = 10;
    mu = 0;
    sig = 0.15;
    guess = [mu sig];

    % Call fmins
    pars = lsqnonlin(@(g) fitcumgauss(g,cons,respc,conN),guess);

    
else
    disp('Cannot continue: Contrast and Coherence were correlated');
    keyboard
end


function SSE = fitcumgauss(guess,x,y,n)
% a = guess(1);
mu = guess(1);
sigma = guess(2);

Est = normcdf(x,mu,sigma);
% Est = a/(sigma*sqrt(2*pi)) * exp( -( (x-mu).^2 ./ (2.*sigma).^2 ) );
% Est = cumsum(Est) ./ sum(Est);
SSE = sum( ((y - Est).^2).*sqrt(n) );
% bar(x, y); hold on;
% plot(x, Est, 'r-'); hold off
% drawnow