function f = cc_rightchoice( adata, fit)

if ieNotDefined('f')
    f = figure;
end
%%
f = figure;

cone = max(adata(:,5)); cohe = max(adata(:,7)); cont = 0.05; coht = 0.1;
conbins = -cone:cont:cone; conrange = -cone+cont/2:cont:cone-cont/2;
cohbins = -cohe:coht:cohe; cohrange = -cohe+coht/2:coht:cohe-coht/2;
clist = brewermap(4,'Oranges');

attend = {'Coherence','Contrast'};
attend_ = [1 2];
types = {'Control','Attend','Unattend'};
types_ = [1 -1 -1];
ccatch = [-1 0 1];

fits = {};

clist = brewermap(9,'PuOr');
for ti = 1:length(types)
    for ai = 1:length(attend)
        flip = [2 1];
        if ccatch(ti)==1
            % this is actually the opposite condition, so flip attend
            data = sel(adata,1,types_(ti)*flip(attend_(ai)));
        else
            data = sel(adata,1,types_(ti)*attend_(ai));
        end
        data = sel(data,9,ccatch(ti));
        subplot(3,2,(ti-1)*2+ai), hold on
        title(sprintf('%s: %s, trials: %i',attend{ai},types{ti},size(data,1)));
        conpeds = unique(data(:,10));
        for ci = 1:length(conpeds)
%             dat = sel(data,10,conpeds(ci));
%             cone = max(dat(:,5)-dat(:,4)); cont = cone/5;
%             cr = -cone+cont/2:cont:cone-cont/2;
            cr = [-1:.01:-.01 .01:.01:1];
%             [mu,std] = buildRcurve(dat(:,8),dat(:,5)-dat(:,4),conbins);
%             plot(conrange,mu,'o','MarkerFaceColor',clist(ci,:),'MarkerEdgeColor',[1 1 1]);
%             errbar(conrange,mu,std,'Color',clist(ci,:));
            [cr,fitcon] = fitCurveCon(fit,conpeds(ci),cr,attend_(ai),ccatch(ti));
            plot(cr,fitcon,'-','Color',clist(ci,:));
        end
        cohpeds = unique(data(:,11));
        for ci = 1:length(cohpeds)
%             dat = sel(data,11,cohpeds(ci));
%             cohe = max(dat(:,7)-dat(:,6));
%             cr = -cohe+coht/2:coht:cohe-coht/2;
            cr = [-1:.01:-.01 .01:.01:1];
            [cr,fitcoh] = fitCurveCoh(fit,cohpeds(ci),cr,attend_(ai),ccatch(ti));
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
        drawPublishAxis
    end
end

function [conr, resp] = fitCurveCon(fit,conped,conrange,attend,ccatch)
%%
if ccatch==1
    % catch trial: adjust the slope of the response function (if necessary)
    if fit.params.conmodel==1
        fit.params.conslope = fit.params.conslope*fit.params.conunatt;
    else
        fit.params.conRmax = fit.params.conRmax*fit.params.conunatt;
    end
end
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
    if ccatch==-1
        diffresp = diffresp * fit.params.beta_control_con_conw;
    elseif ccatch==0
        if isfield(fit.params,'alpha_att_con_conw')
            diffresp = diffresp * (fit.params.beta_control_con_conw + fit.params.alpha_att_con_conw);
        else
            diffresp = diffresp * fit.params.beta_att_con_conw;
        end
    elseif ccatch==1
        if isfield(fit.params,'alpha_unatt_coh_conw')
            diffresp = diffresp * (fit.params.beta_control_coh_conw + fit.params.alpha_unatt_coh_conw);
        else
            diffresp = diffresp * fit.params.beta_unatt_coh_conw;
        end
    end
else
    if ccatch==-1
        diffresp = diffresp * fit.params.beta_control_coh_conw;
    elseif ccatch==0
        if isfield(fit.params,'alpha_att_coh_conw')
            diffresp = diffresp * (fit.params.beta_control_coh_conw + fit.params.alpha_att_coh_conw);
        else
            diffresp = diffresp * fit.params.beta_att_coh_conw;
        end
    elseif ccatch==1
        if isfield(fit.params,'alpha_unatt_con_conw')
            diffresp = diffresp * (fit.params.beta_control_con_conw + fit.params.alpha_unatt_con_conw);
        else
            diffresp = diffresp * fit.params.beta_unatt_con_conw;
        end
    end
end
conr = [-fliplr(conr) conr];
diffresp = [-fliplr(diffresp) diffresp];
diffresp = diffresp + fit.params.bias/2;
if fit.params.poissonNoise
    resp = normcdf(diffresp,0,sqrt(diffresp*fit.params.sigma));
else
    resp = normcdf(diffresp,0,fit.params.sigma);
end
% if mod(length(conrange),2)==0
%     resp = [1-fliplr(resp) resp];
% else
%     resp = [1-fliplr(resp) 0.5 resp];
% end

function [cohr,resp] = fitCurveCoh(fit,cohped,cohrange,attend,ccatch)
%%
if ccatch==1
    if fit.params.cohmodel==1
        fit.params.cohslope = fit.params.cohslope * fit.params.cohunatt;
    else
        fit.params.cohRmax = fit.params.cohRmax*fit.params.cohunatt;
    end
end

cohr = cohrange;
% restrict to values >0 up to pedestal + range = 1
cohr = cohr(cohr>=0);
cohr = cohr(cohped+cohr<=1);
baseresp = cohModel(cohped,fit.params);

highresp = cohModel(cohped+cohr,fit.params);
diffresp = highresp-baseresp;
if attend==1
    if ccatch==-1
        diffresp = diffresp * fit.params.beta_control_coh_cohw;
    elseif ccatch==0
        if isfield(fit.params,'alpha_att_coh_cohw')
            diffresp = diffresp * (fit.params.beta_control_coh_cohw + fit.params.alpha_att_coh_cohw);
        else
            diffresp = diffresp * fit.params.beta_att_coh_cohw;
        end
    elseif ccatch==1
        if isfield(fit.params,'alpha_unatt_con_cohw')
            diffresp = diffresp * (fit.params.beta_control_con_cohw + fit.params.alpha_unatt_con_cohw);
        else
            diffresp = diffresp * fit.params.beta_unatt_con_cohw;
        end
    end
else
    if ccatch==-1
        diffresp = diffresp * fit.params.beta_control_con_cohw;
    elseif ccatch==0
        if isfield(fit.params,'alpha_att_con_cohw')
            diffresp = diffresp * (fit.params.beta_control_con_cohw + fit.params.alpha_att_con_cohw);
        else
            diffresp = diffresp * fit.params.beta_att_con_cohw;
        end
    elseif ccatch==1
        if isfield(fit.params,'alpha_unatt_coh_cohw')
            diffresp = diffresp * (fit.params.beta_control_coh_cohw + fit.params.alpha_unatt_coh_cohw);
        else
            diffresp = diffresp * fit.params.beta_unatt_coh_cohw;
        end
    end
end
% reverse the diffresp to cover negative space
cohr = [-fliplr(cohr) cohr];
diffresp = [-fliplr(diffresp) diffresp];
diffresp = diffresp + fit.params.bias/2;
if fit.params.poissonNoise
    resp = normcdf(diffresp,0,sqrt(diffresp*fit.params.sigma));
else
    resp = normcdf(diffresp,0,fit.params.sigma);
end
% if mod(length(cohrange),2)==0
%     resp = [1-fliplr(resp) resp];
% else
%     resp = [1-fliplr(resp) 0.5 resp];
% end

function out = conModel(con,params)

if params.conmodel==1
    out = params.conslope .* con;
elseif params.conmodel==2
    out = params.conRmax .* ((con.^params.conn) ./ (con.^params.conn + params.conc50.^params.conn));
end

function out = cohModel(coh,params)

if params.cohmodel==1
    out = params.cohslope .* coh;
elseif params.cohmodel==2
    out = params.cohRmax .* ((coh.^params.cohn) ./ (coh.^params.cohn + params.cohc50.^params.cohn));
end