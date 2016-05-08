function allgauss = cc_allcumgauss( adata )
%UNTITLED8 Fit the 6 cumgauss functions for all the conditions:
%
% Attend Contrast:
%   nocatch
%   main
%   catch
%
% Attend Coherence
%   nocatch
%   main
%   catch
%

% NOCATCH
nocatch.con.data = sel(adata,1,2);
[~,fit] = cc_fitCumGauss(nocatch.con.data);
allgauss.nocatch.con.concum = fit.fit.concum;
allgauss.nocatch.con.cohcum = fit.fit.cohcum;
% MAIN
main.con.data = sel(adata,1,-2);
main.con.data = sel(main.con.data,5,0);
% CATCH
ccatch.con.data = sel(adata,1,-2);
ccatch.con.data = sel(ccatch.con.data,5,1);

nocatch.coh.data = sel(adata,1,1);
main.coh.data = sel(adata,1,-1);
main.coh.data = sel(main.coh.data,5,0);
ccatch.coh.data = sel(adata,1,-1);
ccatch.coh.data = sel(ccatch.coh.data,5,1);


