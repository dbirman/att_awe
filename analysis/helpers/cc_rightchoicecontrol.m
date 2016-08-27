function f = cc_rightchoicecontrol( adata, fit)
% Identical to cc_rightchoice, but only uses control trials

adata = adata(adata(:,9)==-1,:);

%%
f = figure;

cone = max(adata(:,5)); cohe = max(adata(:,7)); cont = 0.05; coht = 0.1;
conbins = -cone:cont:cone; conrange = -cone+cont/2:cont:cone-cont/2;
cohbins = -cohe:coht:cohe; cohrange = -cohe+coht/2:coht:cohe-coht/2;
clist = brewermap(4,'Oranges');

attend = {'Coherence','Contrast'};
attend_ = [1 2];

clist = brewermap(9,'PuOr');
for ai = 1:length(attend)
    data = sel(adata,1,attend_(ai));
    subplot(1,2,ai), hold on
    title(sprintf('%s, n=%i',attend{ai},size(data,1)));
    conpeds = unique(data(:,10));
    for ci = 1:length(conpeds)
        %             dat = sel(data,10,conpeds(ci));
        %             cone = max(dat(:,5)-dat(:,4)); cont = cone/5;
        %             cr = -cone+cont/2:cont:cone-cont/2;
        cr = [-1:.01:-.01 .01:.01:1];
        %             [mu,std] = buildRcurve(dat(:,8),dat(:,5)-dat(:,4),conbins);
        %             plot(conrange,mu,'o','MarkerFaceColor',clist(ci,:),'MarkerEdgeColor',[1 1 1]);
        %             errbar(conrange,mu,std,'Color',clist(ci,:));
        [cr,fitcon] = fitCurveCon(fit,conpeds(ci),cr,attend_(ai));
        plot(cr,fitcon,'-','Color',clist(ci,:));
    end
    cohpeds = unique(data(:,11));
    for ci = 1:length(cohpeds)
        %             dat = sel(data,11,cohpeds(ci));
        %             cohe = max(dat(:,7)-dat(:,6));
        %             cr = -cohe+coht/2:coht:cohe-coht/2;
        cr = [-1:.01:-.01 .01:.01:1];
        [cr,fitcoh] = fitCurveCoh(fit,cohpeds(ci),cr,attend_(ai));
        plot(cr,fitcoh,'-','Color',clist(10-ci,:));
    end
    [conmu,constd] = buildRcurve(data(:,8),data(:,5)-data(:,4),conbins);
    [cohmu,cohstd] = buildRcurve(data(:,8),data(:,7)-data(:,6),cohbins);
    plot(conrange,conmu,'o','MarkerFaceColor',clist(1,:),'MarkerEdgeColor',[1 1 1]);
    errbar(conrange,conmu,constd,'Color',clist(1,:));
    plot(cohrange,cohmu,'o','MarkerFaceColor',clist(9,:),'MarkerEdgeColor',[1 1 1]);
    errbar(cohrange,cohmu,cohstd,'Color',clist(9,:));
    %         try
    %         [~,fit] = cc_fitCumGauss(data);
    %         fits{end+1} = fit;
    %         plot(fit.fit.x,fit.fit.concum,'-','Color',clist(1,:));
    %         plot(fit.fit.x,fit.fit.cohcum,'-','Color',clist(3,:));
    %         title(sprintf('Task: %s Cond: %s Alpha: %0.2f',attend{ai},types{ti},fit.params(5)))
    %         catch
    %         end
    axis([-1 1 0 1]);
    
    % set(ax,'XTickLabel',{'-3pi','-2pi','-pi','0','pi','2pi','3pi'})
    set(gca,'XTickLabel',{'-100','-50','0','50','100'});
    set(gca,'YTick',[0 0.5 1],'YTickLabel',{'0','','100%'});
    xlabel('Stimulus Strength (Right - Left, %)');
    ylabel('Right Choice Probability');
    drawPublishAxis
end

function [conr, resp] = fitCurveCon(fit,conped,conrange,attend)

% get the contrast range and compute the model response
% conr = conrange(conrange>=0); % only model additive changes, flip for the negative
conr = conrange;
conr = conr(conr>=0);
conr = conr(conped+conr<1);
baseresp = conModel(conped,fit.params);
highresp = conModel(conped+conr,fit.params);
diffresp = highresp-baseresp;
% adjust the response by the condition, so using alphacon
if attend==2
    % CONTRAST
    diffresp = diffresp * fit.params.beta_control_con_conw;
else
    diffresp = diffresp * fit.params.beta_control_coh_conw;
end
conr = [-fliplr(conr) conr];
diffresp = [-fliplr(diffresp) diffresp];
diffresp = diffresp + fit.params.bias;
if fit.params.poissonNoise
    resp = normcdf(diffresp,0,sqrt(abs(diffresp*fit.params.sigma)));
else
    resp = normcdf(diffresp,0,fit.params.sigma);
end

function [cohr,resp] = fitCurveCoh(fit,cohped,cohrange,attend)
%%

cohr = cohrange;
% restrict to values >0 up to pedestal + range = 1
cohr = cohr(cohr>=0);
cohr = cohr(cohped+cohr<=1);
baseresp = cohModel(cohped,fit.params);

highresp = cohModel(cohped+cohr,fit.params);
diffresp = highresp-baseresp;
if attend==1
    diffresp = diffresp * fit.params.beta_control_coh_cohw;
else
    diffresp = diffresp * fit.params.beta_control_con_cohw;
end
% reverse the diffresp to cover negative space
cohr = [-fliplr(cohr) cohr];
diffresp = [-fliplr(diffresp) diffresp];
diffresp = diffresp + fit.params.bias/2;
if fit.params.poissonNoise
    resp = normcdf(diffresp,0,sqrt(abs(diffresp*fit.params.sigma)));
else
    resp = normcdf(diffresp,0,fit.params.sigma);
end
% if mod(length(cohrange),2)==0
%     resp = [1-fliplr(resp) resp];
% else
%     resp = [1-fliplr(resp) 0.5 resp];
% end
