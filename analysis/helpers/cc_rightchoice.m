function cc_rightchoice( adata, f )

if ieNotDefined('f')
    f = figure;
end
%%
figure(f);
clf
cone = max(adata(:,2)); cohe = max(adata(:,3)); cont = 0.05; coht = 0.1;
conbins = -cone:cont:cone; conrange = -cone+cont/2:cont:cone-cont/2;
cohbins = -cohe:coht:cohe; cohrange = -cohe+coht/2:coht:cohe-coht/2;
clist = brewermap(4,'Oranges');

attend = {'Coherence','Contrast'};
attend_ = [1 2];
types = {'Control','Attend','Unattend'};
types_ = [1 -1 -1];
ccatch = [-1 0 1];

fits = {};

clist = brewermap(3,'PuOr');
for ti = 1:length(types)
    for ai = 1:length(attend)
        data = sel(adata,1,types_(ti)*attend_(ai));
        data = sel(data,5,ccatch(ti));
        subplot(3,2,(ti-1)*2+ai), hold on
        [conmu,constd] = buildRcurve(data(:,4),data(:,2),conbins);
        [cohmu,cohstd] = buildRcurve(data(:,4),data(:,3),cohbins);
        plot(conrange,conmu,'o','MarkerFaceColor',clist(1,:),'MarkerEdgeColor',[1 1 1]);
        errbar(conrange,conmu,constd,'Color',clist(1,:));
        plot(cohrange,cohmu,'o','MarkerFaceColor',clist(3,:),'MarkerEdgeColor',[1 1 1]);
        errbar(cohrange,cohmu,cohstd,'Color',clist(3,:));
%         try
        [~,fit] = cc_fitCumGauss(data);
        fits{end+1} = fit;
        plot(fit.fit.x,fit.fit.concum,'-','Color',clist(1,:));
        plot(fit.fit.x,fit.fit.cohcum,'-','Color',clist(3,:));
        title(sprintf('Task: %s Cond: %s Alpha: %0.2f',attend{ai},types{ti},fit.params(5)))
%         catch
%         end
        axis([-1 1 0 1]);
    end
end

