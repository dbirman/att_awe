function cc_rightChoicePlots(dofit)


clist2 = brewermap(3,'Oranges');
clist1 = brewermap(3,'Blues');

% allData = loadAllData(sid);

fullData = loadFullData;

acon = fullData.data(:,5)'-fullData.data(:,4)';
acoh = fullData.data(:,7)'-fullData.data(:,6)';
aresp = fullData.data(:,12)';
atask = fullData.data(:,3)';
acatch = fullData.data(:,11)';
anocatch = fullData.data(:,13)';

tasks = {'coherence','contrast'};
catchs = {'catch','nocatch'};
cueds = {'control','cued','miscued'};

name = {};
range = {};
values = {};

for si = 1:length(fullData.subjects)
    sid = fullData.subjects{si};
    
    subj = (fullData.data(:,14)==mrStr2num(strrep(sid,'s','')))';
    
    disp(sprintf('Subject %s completed %i trials',sid,sum(subj)));
    figure;
    for nci = 0:1
        nocatchs = catchs{nci+1};
        
        for ti = 1:length(tasks)
            ctask = tasks{ti};
            
            % pick trials
            choose = logical(subj .* (anocatch==nci) .* (atask==ti));
            % select from atrials
            con = acon(choose); coh = acoh(choose); resp = aresp(choose)-1;
            task = atask(choose); ccatch = logical(acatch(choose));
            
            cone = .3; cohe = .8; cont = 0.025; coht = 0.1;
            conbins = -cone:cont:cone; conrange = -cone+cont/2:cont:cone-cont/2;
            cohbins = -cohe:coht:cohe; cohrange = -cohe+coht/2:coht:cohe-coht/2;
            % if this is a catch trial run we should separate out the effect of
            % cued and miscued trials separately
            if nci==0

                if dofit==1
                    X = [coh',con',(coh.*ccatch)'];
                    Y = resp';
                    [B,dev,stats] = glmfit(X,Y,'binomial','link','logit');
                elseif dofit==2
                    % cum gaussian
                    con1 = con(ccatch); coh1 = coh(ccatch); resp1 = resp(ccatch);
                    con2 = con(~ccatch); coh2 = coh(~ccatch); resp2 = resp(~ccatch);
                    pars_miscued = cc_fitCumGauss(con1,coh1,resp1,0);
                    pars_cued = cc_fitCumGauss(con2,coh2,resp2,0);
                end

                % get actual data points for CUED trials
                conBinned = binData(resp(ccatch==0),con(ccatch==0),conbins);
                conMeanCued = cellfun(@nanmean,conBinned);
                conStdCued = cellfun(@(A) 0,conBinned);
                % coherence data points
                cohBinned = binData(resp(ccatch==0),coh(ccatch==0),cohbins);
                cohMeanCued = cellfun(@nanmean,cohBinned);
                cohStdCued = cellfun(@(A) 0,cohBinned);

                % get actual data points for MISCUED trials
                conBinned2 = binData(resp(ccatch==1),con(ccatch==1),conbins);
                conMeanMiscued = cellfun(@nanmean,conBinned2);
                conStdMiscued = cellfun(@(A) 0,conBinned2);
                % coherence data points
                cohBinned2 = binData(resp(ccatch==1),coh(ccatch==1),cohbins);
                cohMeanMiscued = cellfun(@nanmean,cohBinned2);
                cohStdMiscued = cellfun(@(A) 0,cohBinned2);

                subplot(3,2,ti);
                hold on

                % plot cued trials

                % P(R) = Icoh + Icon + iscatch
                if dofit==1 
                    y = glmval(B,[zeros(61,1),(-.3:.01:.3)',zeros(61,1)],'logit');
                elseif dofit==2
                    y = normcdf(-.3:.01:.3,pars_cued(1),pars_cued(2));
                end
                if dofit>0
                    plot(-.3:.01:.3,y,'-','Color',clist2(3,:));
                    name{end+1} = 'contrast cued (catch)';
                    range{end+1} = -.3:.01:.3;
                    values{end+1} = y;
                end
                %                 errorbar(conrange,conMeanCued,conStdCued,'or');
                plot(conrange,conMeanCued,'o','Color',clist2(3,:));
%                 plot(conrange,conMeanCued,'-r');

                if dofit==1
                    y = glmval(B,[(-.85:.01:.85)',zeros(171,1),zeros(171,1)],'logit');
                elseif dofit==2
                    y = normcdf(-.85:.01:.85,pars_cued(3),pars_cued(4));
                end
                if dofit>0
                    plot(-.85:.01:.85,y,'-','Color',clist1(3,:))
                    name{end+1} = 'coherence cued (catch)';
                    range{end+1} = -.85:.01:.85;
                    values{end+1} = y; 
                end
%                 errorbar(cohrange,cohMeanCued,cohStdCued,'o','Color',clist1(3,:));
                plot(cohrange,cohMeanCued,'o','Color',clist1(3,:));
%                 plot(cohrange,cohMeanCued,'-','Color',clist1(3,:));

                title(sprintf('Cued: Response for %s, %s',ctask,nocatchs));

                % plot miscued trials
                subplot(3,2,4+ti);
                hold on
                if dofit==1
                    y = glmval(B,[zeros(61,1),(-.3:.01:.3)',zeros(61,1)],'logit');
                elseif dofit==2
                    y = normcdf(-.3:.01:.3,pars_miscued(1),pars_miscued(2));
                end
                if dofit>0
                    plot(-.3:.01:.3,y,'-r');
                end
%                 errorbar(conrange,conMeanMiscued,conStdMiscued,'or');
                plot(conrange,conMeanMiscued,'o','Color',clist2(3,:));
%                 plot(conrange,conMeanMiscued,'-r');

                if dofit==1
                    y = glmval(B,[(-.85:.01:.85)',zeros(171,1),(-.85:.01:.85)'],'logit');
                elseif dofit==2
                    y = normcdf(-.85:.01:.85,pars_miscued(3),pars_miscued(4));
                end
                if dofit>0
                    plot(-.85:.01:.85,y,'-','Color',clist1(3,:))
                end
%                 errorbar(cohrange,cohMeanMiscued,cohStdMiscued,'o','Color',clist1(3,:));
                plot(cohrange,cohMeanMiscued,'o','Color',clist1(3,:));
%                 plot(cohrange,cohMeanMiscued,'-','Color',clist1(3,:));
                title(sprintf('Mis-Cued: Response for %s, %s',ctask,nocatchs));
            else
                % run logistic regression
                if dofit==1
                    X = [coh',con'];
                    Y = resp';
                    [B,dev,stats] = glmfit(X,Y,'binomial','link','logit');
                elseif dofit==2
                    pars = cc_fitCumGauss(con,coh,resp,0);
                end

                % get actual data points
                conBinned = binData(resp,con,conbins);
                conMean = cellfun(@nanmean,conBinned);
                conLength = cellfun(@length,conBinned);
%                 conStd = cellfun(@(A) 0,conBinned);
                conMean = conMean(conLength>10);
%                 conStd = conStd(conLength>5);
                % coherence data points
                cohBinned = binData(resp,coh,cohbins);
                cohMean = cellfun(@nanmean,cohBinned);
                cohLength = cellfun(@length,cohBinned);
%                 cohStd = cellfun(@(A) 0,cohBinned);
                cohMean = cohMean(cohLength>10); %cohStd = cohStd(cohLength>5);

                subplot(3,2,2+ti);
                hold on
                if dofit==1
                    y = glmval(B,[zeros(61,1),(-.3:.01:.3)'],'logit');
                elseif dofit==2
                    y = normcdf(-.3:.01:.3,pars(1),pars(2));
                end  
                if dofit>0
                    plot(-.3:.01:.3,y,'-','Color',clist2(3,:));  
                end
%                 errorbar(conrange,conMean,conStd,'or');
%                 plot(conrange(conLength>10),conMean,'-r');
                plot(conrange(conLength>10),conMean,'or');
                if dofit==1
                    y = glmval(B,[(-.85:.01:.85)',zeros(171,1)],'logit');
                elseif dofit==2
                    y = normcdf(-.85:.01:.85,pars(3),pars(4));
                end
                if dofit>0
                    plot(-.85:.01:.85,y,'-','Color',clist1(3,:));
                end
%                 errorbar(cohrange,cohMean,cohStd,'o','Color',clist1(3,:));
%                 plot(cohrange(cohLength>10),cohMean,'-','Color',clist1(3,:));
                plot(cohrange(cohLength>10),cohMean,'o','Color',clist1(3,:));
                axis([-.85 .85 0 1]);
    %             legend({'Contrast','Coherence'});
                title(sprintf('Control: Response for %s, %s',ctask,nocatchs));
            end
        end
    end
    dir = fullfile('~/proj/att_awe/analysis/figures',sid);
    if ~isdir(dir), mkdir(dir); end
    fname = fullfile(dir,'rightChoice.pdf');
    print(fname,'-dpdf');
end