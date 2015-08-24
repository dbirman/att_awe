function cc_rightChoicePlots(sid)

allData = loadAllData(sid);

I = allData.behav.I;

tasks = {'coherence','contrast'};
catchs = {'catch','nocatch'};
cueds = {'control','cued','miscued'};

        figure
for ci = 1:length(catchs)
    catchc = catchs{ci};
    for ti = 1:length(tasks)
        ctask = tasks{ti};
        
        con = I.(catchc).con;
        coh = I.(catchc).coh;
        side = I.(catchc).side;
        task = I.(catchc).task;
        
        conbins = -.3:.05:.3; conrange = -.275:.05:.275;
        cohbins = -.8:.1:.8; cohrange = -.75:.1:.75;
        % if this is a catch trial run we should separate out the effect of
        % cued and miscued trials separately
        if strcmp(catchc,'catch')
            ccatch = I.(catchc).catch;
            conv = con(task==ti); cohv = coh(task==ti); catv = ccatch(task==ti);
            X = [cohv',conv',(cohv.*catv)'];
            Y = side(task==ti)'-1;
            [B,dev,stats] = glmfit(X,Y,'binomial','link','logit');
            
            % one more beta value here, we have to deal with that
            % get actual data points
            conBinned = binData(side(task==ti)-1,con(task==ti),conbins);
            conMean = cellfun(@mean,conBinned);
            conStd = cellfun(@(A) std(A)/sqrt(length(A))*1.96,conBinned);
            % coherence data points
            cohBinned = binData(side(task==ti)-1,coh(task==ti),cohbins);
            cohMean = cellfun(@mean,cohBinned);
            cohStd = cellfun(@(A) std(A)/sqrt(length(A))*1.96,cohBinned);
            
            subplot(3,2,(ci-1)*2+ti);
            hold on
            % plot cued trials
            
            % plot main trials (on catch runs)
            % P(R) = Icoh + Icon + iscatch
            y = glmval(B,[zeros(61,1),(-.3:.01:.3)',zeros(61,1)],'logit');
            plot(-.3:.01:.3,y,'-r');
            
            y = glmval(B,[(-.85:.01:.85)',zeros(171,1),zeros(171,1)],'logit');
            plot(-.85:.01:.85,y,'-b')
            
            title(sprintf('Cued: Response for %s, %s',ctask,catchc));
            
            subplot(3,2,4+ti);
            hold on
            y = glmval(B,[zeros(61,1),(-.3:.01:.3)',zeros(61,1)],'logit');
            plot(-.3:.01:.3,y,'-r');
            
            y = glmval(B,[(-.85:.01:.85)',zeros(171,1),(-.85:.01:.85)'],'logit');
            plot(-.85:.01:.85,y,'-b')
            title(sprintf('Mis-Cued: Response for %s, %s',ctask,catchc));
            % plot miscued trials
        else
            % run logistic regression
            X = [coh(task==ti)',con(task==ti)'];
            Y = side(task==ti)'-1;
            [B,dev,stats] = glmfit(X,Y,'binomial','link','logit');
        
            % get actual data points
            conBinned = binData(side(task==ti)-1,con(task==ti),conbins);
            conMean = cellfun(@mean,conBinned);
            conStd = cellfun(@(A) std(A)/sqrt(length(A))*1.96,conBinned);
            % coherence data points
            cohBinned = binData(side(task==ti)-1,coh(task==ti),cohbins);
            cohMean = cellfun(@mean,cohBinned);
            cohStd = cellfun(@(A) std(A)/sqrt(length(A))*1.96,cohBinned);

            subplot(3,2,(ci-1)*2+ti);
            hold on
            y = glmval(B,[zeros(61,1),(-.3:.01:.3)'],'logit');
            plot(-.3:.01:.3,y,'-r')    
            errorbar(conrange,conMean,conStd,'-r');
            y = glmval(B,[(-.85:.01:.85)',zeros(171,1)],'logit');
            plot(-.85:.01:.85,y,'-b')
            errorbar(cohrange,cohMean,cohStd,'-b');
            axis([-.85 .85 0 1]);
            legend({'Contrast','Coherence'});
            title(sprintf('Control: Response for %s, %s',ctask,catchc));
        end
    end
end