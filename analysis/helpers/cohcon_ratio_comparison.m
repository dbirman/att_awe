
%% Compare responses
x = 0:.01:1;
% 
con = zeros(length(nSIDs),4,101);
coh = zeros(length(nSIDs),4,101);
for si = 1:length(nSIDs)
    for di = 1:4
        mfit = sfits{si}{di};
        con(si,di,:) = conModel(x,mfit.roifit{1}.params);
        coh(si,di,:) = cohModel(x,mfit.roifit{8}.params);
    end
end
    
%% plot
rat = con./coh;
rat = squeeze(mean(rat));

con_ = squeeze(mean(con));
coh_ = squeeze(mean(coh));

figure; hold on
cmap = brewermap(7,'Oranges');
cmapp = brewermap(7,'Purples');

for di = 1:4
    p(di) = plot(x,con_(di,:),'Color',cmap(di,:));
    plot(x,coh_(di,:),'Color',cmapp(di,:));
end
legend(p,dataopts);

figure; hold on
plot(x,rat);
axis([0 1 0 25]);
legend(dataopts);