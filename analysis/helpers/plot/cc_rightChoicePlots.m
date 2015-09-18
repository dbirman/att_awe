function cc_rightChoicePlots(dofit)


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
            con = acon(choose); coh = acoh(choose); side = aresp(choose)-1;
            task = atask(choose); ccatch = acatch(choose);
            
            cone = .3; cohe = .8; cont = 0.025; coht = 0.1;
            conbins = -cone:cont:cone; conrange = -cone+cont/2:cont:cone-cont/2;
            cohbins = -cohe:coht:cohe; cohrange = -cohe+coht/2:coht:cohe-coht/2;
            % if this is a catch trial run we should separate out the effect of
            % cued and miscued trials separately
            if nci==0

                if dofit
                X = [cohv',conv',(cohv.*catv)'];
                Y = side(task==ti)'-1;
                [B,dev,stats] = glmfit(X,Y,'binomial','link','logit');
                end

                % get actual data points for CUED trials
                conBinned = binData(side(ccatch==0),con(ccatch==0),conbins);
                conMeanCued = cellfun(@nanmean,conBinned);
                conStdCued = cellfun(@(A) 0,conBinned);
                % coherence data points
                cohBinned = binData(side(ccatch==0),coh(ccatch==0),cohbins);
                cohMeanCued = cellfun(@nanmean,cohBinned);
                cohStdCued = cellfun(@(A) 0,cohBinned);

                % get actual data points for MISCUED trials
                conBinned2 = binData(side(ccatch==1),con(ccatch==1),conbins);
                conMeanMiscued = cellfun(@nanmean,conBinned2);
                conStdMiscued = cellfun(@(A) 0,conBinned2);
                % coherence data points
                cohBinned2 = binData(side(ccatch==1),coh(ccatch==1),cohbins);
                cohMeanMiscued = cellfun(@nanmean,cohBinned2);
                cohStdMiscued = cellfun(@(A) 0,cohBinned2);

                subplot(3,2,ti);
                hold on

                % plot cued trials

                % P(R) = Icoh + Icon + iscatch
                if dofit 
                    y = glmval(B,[zeros(61,1),(-.3:.01:.3)',zeros(61,1)],'logit');
                    plot(-.3:.01:.3,y,'-r');
                end
%                 errorbar(conrange,conMeanCued,conStdCued,'*r');
                plot(conrange,conMeanCued,'*r');
                plot(conrange,conMeanCued,'-r');

                if dofit 
                    y = glmval(B,[(-.85:.01:.85)',zeros(171,1),zeros(171,1)],'logit');
                    plot(-.85:.01:.85,y,'-b')
                end
%                 errorbar(cohrange,cohMeanCued,cohStdCued,'*b');
                plot(cohrange,cohMeanCued,'*b');
                plot(cohrange,cohMeanCued,'-b');

                title(sprintf('Cued: Response for %s, %s',ctask,nocatchs));

                % plot miscued trials
                subplot(3,2,4+ti);
                hold on
                if dofit 
                    y = glmval(B,[zeros(61,1),(-.3:.01:.3)',zeros(61,1)],'logit');
                    plot(-.3:.01:.3,y,'-r');
                end
%                 errorbar(conrange,conMeanMiscued,conStdMiscued,'*r');
                plot(conrange,conMeanMiscued,'*r');
                plot(conrange,conMeanMiscued,'-r');

                if dofit 
                    y = glmval(B,[(-.85:.01:.85)',zeros(171,1),(-.85:.01:.85)'],'logit');
                    plot(-.85:.01:.85,y,'-b')
                end
%                 errorbar(cohrange,cohMeanMiscued,cohStdMiscued,'*b');
                plot(cohrange,cohMeanMiscued,'*b');
                plot(cohrange,cohMeanMiscued,'-b');
                title(sprintf('Mis-Cued: Response for %s, %s',ctask,nocatchs));
            else
                % run logistic regression
                if dofit
                    X = [coh(task==ti)',con(task==ti)'];
                    Y = side(task==ti)'-1;
                    [B,dev,stats] = glmfit(X,Y,'binomial','link','logit');
                end

                % get actual data points
                conBinned = binData(side,con,conbins);
                conMean = cellfun(@nanmean,conBinned);
                conLength = cellfun(@length,conBinned);
%                 conStd = cellfun(@(A) 0,conBinned);
                conMean = conMean(conLength>10);
%                 conStd = conStd(conLength>5);
                % coherence data points
                cohBinned = binData(side,coh,cohbins);
                cohMean = cellfun(@nanmean,cohBinned);
                cohLength = cellfun(@length,cohBinned);
%                 cohStd = cellfun(@(A) 0,cohBinned);
                cohMean = cohMean(cohLength>10); %cohStd = cohStd(cohLength>5);

                subplot(3,2,2+ti);
                hold on
                if dofit 
                    y = glmval(B,[zeros(61,1),(-.3:.01:.3)'],'logit');
                    plot(-.3:.01:.3,y,'-r')  
                end  
%                 errorbar(conrange,conMean,conStd,'*r');
                plot(conrange(conLength>10),conMean,'-r');
                plot(conrange(conLength>10),conMean,'*r');
                if dofit 
                    y = glmval(B,[(-.85:.01:.85)',zeros(171,1)],'logit');
                    plot(-.85:.01:.85,y,'-b')
                end
%                 errorbar(cohrange,cohMean,cohStd,'*b');
                plot(cohrange(cohLength>10),cohMean,'-b');
                plot(cohrange(cohLength>10),cohMean,'*b');
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