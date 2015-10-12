function cc_neuralFigures(neural,name,sid,df, nf)

dir = fullfile('~/proj/att_awe/analysis/figures',sid);
dir1 = fullfile('~/proj/att_awe/analysis/figures',sid,'single_panels');
if ~isdir(dir), mkdir(dir); end
if ~isdir(dir1), mkdir(dir1); end
if ~isdir(fullfile(dir1,'deconv')), mkdir(fullfile(dir1,'deconv')); end
if ~isdir(fullfile(dir1,'naka')), mkdir(fullfile(dir1,'naka')); end
if ~isdir(fullfile(dir1,'roi')), mkdir(fullfile(dir1,'roi')); end

roi_effects = []; % roi x feature x task x effect

rois = neural.shortROIs;
deconvo = neural.(name).deconvo;
fits = neural.(name).fits;
% Contrast Response at each Coherence
tasks = {'Cued Coherence','Cued Contrast'};
%% Deconv Figures
if df
for ri = 1:length(rois)
    [conVals, cohVals, cuedTask] = parseNames(neural.SCM.(name).(rois{ri}).stimNames); 
    ucoh = unique(cohVals);
    ucon = unique(conVals);
    clist2 = brewermap(length(ucoh),'Oranges');
    clist1 = brewermap(length(ucon),'Blues');
    if ~length(ucoh)==length(ucon)
        disp('This code will break.. fix it up');keyboard
    end
    roi = rois{ri};
    dfs(1) = figure('name',sprintf('%s: %s',tasks{1},roi)); dfs(2) = figure('name',sprintf('%s: %s',tasks{2},roi));
%     fit = fits{ri}.fit;
    fit = deconvo{ri};
    
    whitebg(dfs(1),'k')
    whitebg(dfs(2),'k')
    
    for task = 1:2
        title(sprintf('Task %s',tasks{task}));
        for ui = 1:length(ucoh)
            ccoh = ucoh(ui);
            figure(dfs(task));
            subplot(length(ucoh),2,2*ui-1);
            hold on
            % draw change in contrast at each coherence
            for ci = 1:length(ucon)
                ccon = ucon(ci);
                % get the ehdr for con * coh * task
                pcon = conVals==ccon;
                pcoh = cohVals==ccoh;
                ptask = cuedTask==task;
                ehdr = fit.ehdr(logical(pcon.*pcoh.*ptask),:);
                se = fit.ehdrste(logical(pcon.*pcoh.*ptask),:);
                if ~isempty(ehdr)
                plot(fit.time,ehdr,'Color',clist1(ci,:));
                end
%                 errorbar(ehdr,se,'Color',clist1(ci,:));
            end
            title(sprintf('Coherence: %0.02f',ccoh));
            [f,c] = sp_figure; figure(f), hold on
            
            % REPEAT EVERYTHING FOR SP
            
            % draw change in contrast at each coherence
            for ci = 1:length(ucon)
                ccon = ucon(ci);
                % get the ehdr for con * coh * task
                pcon = conVals==ccon;
                pcoh = cohVals==ccoh;
                ptask = cuedTask==task;
                ehdr = fit.ehdr(logical(pcon.*pcoh.*ptask),:);
                se = fit.ehdrste(logical(pcon.*pcoh.*ptask),:);
                if ~isempty(ehdr)
                plot(fit.time,ehdr,'Color',clist1(ci,:));
                end
%                 errorbar(ehdr,se,'Color',clist1(ci,:));
            end
            title(sprintf('Coherence: %0.02f',ccoh));
            a = axis();
            axis([min(fit.time) max(fit.time) a(3) a(4)]);
            fname = fullfile(dir1,'deconv',sprintf('Deconv_%s_Coh_%i',roi,ccoh*100));
            print(fname,'-dpdf');
            close(f)
        end
        for ui = 1:length(ucon)
            ccon = ucon(ui);
            figure(dfs(task));
            subplot(length(ucon),2,2*ui);
            hold on
            % draw change in contrast at each coherence
            for ci = 1:length(ucoh)
                ccoh = ucoh(ci);
                % get the ehdr for con * coh * task
                pcon = conVals==ccon;
                pcoh = cohVals==ccoh;
                ptask = cuedTask==task;
                ehdr = fit.ehdr(logical(pcon.*pcoh.*ptask),:);
                se = fit.ehdrste(logical(pcon.*pcoh.*ptask),:);
                if ~isempty(ehdr)
                plot(fit.time,ehdr,'Color',clist2(ci,:));
                end
%                 errorbar(ehdr,se,'Color',clist2(ci,:));
            end
            title(sprintf('Contrast: %0.02f',ccon));
                        
            [f,c] = sp_figure; figure(f), hold on
            % draw change in contrast at each coherence
            for ci = 1:length(ucoh)
                ccoh = ucoh(ci);
                % get the ehdr for con * coh * task
                pcon = conVals==ccon;
                pcoh = cohVals==ccoh;
                ptask = cuedTask==task;
                ehdr = fit.ehdr(logical(pcon.*pcoh.*ptask),:);
                se = fit.ehdrste(logical(pcon.*pcoh.*ptask),:);
                if ~isempty(ehdr)
                plot(fit.time,ehdr,'Color',clist2(ci,:));
                end
%                 errorbar(ehdr,se,'Color',clist2(ci,:));
            end
            title(sprintf('Contrast: %0.02f',ccon));
            axis([min(fit.time) max(fit.time) a(3) a(4)]);
            fname = fullfile(dir1,'deconv',sprintf('Deconv_%s_Con_%i',roi,ccon*100));
            print(fname,'-dpdf');
            close(f)

        end
        figure(dfs(task))
        axis([min(fit.time) max(fit.time) a(3) a(4)]);
        fname = fullfile(dir,sprintf('Deconv_%s',roi));
        print(fname,'-dpdf');
        close(gcf);

    end
end
end

%% Amplitude Figures
if nf
f1 = figure;
f2 = figure;
for ri = 1:length(rois)
    [conVals, cohVals, cuedTask] = parseNames(neural.SCM.(name).(rois{ri}).stimNames); 
    ucoh = unique(cohVals);
    ucon = unique(conVals);
    clist2 = brewermap(length(ucoh),'Oranges');
    clist1 = brewermap(length(ucon),'Blues');
    if ~length(ucoh)==length(ucon)
        disp('This code will break.. fix it up');keyboard
    end
    roi = rois{ri};     
    fit = fits{ri};
    N = [];
    for si = 1:length(neural.SCM.(name).(roi).stimVol)
        N(end+1) = length(neural.SCM.(name).(roi).stimVol{si});
    end
    conVals = conVals(cuedTask<=2);
    cohVals = cohVals(cuedTask<=2);
    amps = fit.amplitude(cuedTask<=2);
    ampse = fit.amplitudeSTE(cuedTask<=2);
    N = N(cuedTask<=2);
    cuedTask = cuedTask(cuedTask<=2);

    ucon = unique(conVals);
    ucoh = unique(cohVals);

    mat{1} = zeros(length(ucon),length(ucoh));
    mat{2} = zeros(length(ucon),length(ucoh));
    matse{1} = zeros(length(ucon),length(ucoh));
    matse{2} = zeros(length(ucon),length(ucoh));
    Nm{1} = zeros(length(ucon),length(ucoh));
    Nm{2} = zeros(length(ucon),length(ucoh));

    for ampi = 1:length(amps)
        mat{cuedTask(ampi)}(find(ucon==conVals(ampi)),find(ucoh==cohVals(ampi))) = amps(ampi);
        matse{cuedTask(ampi)}(find(ucon==conVals(ampi)),find(ucoh==cohVals(ampi))) = ampse(ampi);
        Nm{cuedTask(ampi)}(find(ucon==conVals(ampi)),find(ucoh==cohVals(ampi))) = N(ampi);
    end

    neural.fit.(name).(roi).mat = mat;
    neural.fit.(name).(roi).matse = matse;
    neural.fit.(name).(roi).Nm = Nm;

    figure(f1)
    subplot(length(rois),2,(ri-1)*2+1); hold on
    % X axis = contrast
    legVals = {};
    for i = 1:4
        plot(ucon,mat{1}(:,i),'-','Color',clist1(i,:));
    end
    for i = 1:4
        plot(ucon,mat{2}(:,i),'-','Color',clist2(i,:));
    end
%     legend(legVals);
    for i = 1:4
        errorbar(ucon,mat{1}(:,i),matse{1}(:,i),'o','Color',clist1(i,:));
        errorbar(ucon,mat{2}(:,i),matse{2}(:,i),'o','Color',clist2(i,:));
    end
    title(sprintf('%s Orange = Attend Con; Blue = Attend Coh',roi));
    xlabel('Contrast');
    ylabel('Response Amplitude');
    miv = min([mat{1}(:) - matse{1}(:);mat{2}(:) - matse{2}(:)])-0.05;
    mav = max([mat{1}(:) + matse{1}(:);mat{2}(:) + matse{2}(:)])+0.05;
    axis([min(ucon)-.05 max(ucon)+.05 miv mav]);
%     set(ax,'YTick',[-1 -0.5 0 0.5 1])

    % REPEAT FOR SP
    [s,f] = sp_figure(); hold on
        for i = 1:4
        plot(ucon,mat{1}(:,i),'-','Color',clist1(i,:));
    end
    for i = 1:4
        plot(ucon,mat{2}(:,i),'-','Color',clist2(i,:));
    end
%     legend(legVals);
    for i = 1:4
        errorbar(ucon,mat{1}(:,i),matse{1}(:,i),'o','Color',clist1(i,:));
        errorbar(ucon,mat{2}(:,i),matse{2}(:,i),'o','Color',clist2(i,:));
    end
    title(sprintf('%s Orange = Attend Con; Blue = Attend Coh',roi));
    xlabel('Contrast');
    ylabel('Response Amplitude');
    miv = min([mat{1}(:) - matse{1}(:);mat{2}(:) - matse{2}(:)])-0.05;
    mav = max([mat{1}(:) + matse{1}(:);mat{2}(:) + matse{2}(:)])+0.05;
    axis([min(ucon)-.05 max(ucon)+.05 miv mav]);
    fname = fullfile(dir1,'naka',sprintf('range_contrast_%s',roi));
    print(fname,'-dpdf');
    close(s)
    
    % X axis = coherence
    subplot(length(rois),2,(ri-1)*2+2), hold on
    for i = 1:4
        plot(ucoh,mat{1}(i,:),'-','Color',clist1(i,:));
        errorbar(ucoh,mat{1}(i,:),matse{1}(i,:),'o','Color',clist1(i,:));
        plot(ucoh,mat{2}(i,:),'-','Color',clist2(i,:));
        errorbar(ucoh,mat{2}(i,:),matse{2}(i,:),'o','Color',clist2(i,:));
    end
    title(sprintf('%s Orange = Attend Con; Blue = Attend Coh',roi));
    xlabel('Coherence');
    ylabel('Response Amplitude');
    axis([min(ucoh)-.05 max(ucoh)+.05 miv mav]);

    % REPEAT FOR SP
    [s,f] = sp_figure(); hold on
    for i = 1:4
        plot(ucoh,mat{1}(i,:),'-','Color',clist1(i,:));
        errorbar(ucoh,mat{1}(i,:),matse{1}(i,:),'o','Color',clist1(i,:));
        plot(ucoh,mat{2}(i,:),'-','Color',clist2(i,:));
        errorbar(ucoh,mat{2}(i,:),matse{2}(i,:),'o','Color',clist2(i,:));
    end
    title(sprintf('%s Orange = Attend Con; Blue = Attend Coh',roi));
    xlabel('Coherence');
    ylabel('Response Amplitude');
    axis([min(ucoh)-.05 max(ucoh)+.05 miv mav]);
    fname = fullfile(dir1,'naka',sprintf('range_coherence_%s',roi));
    print(fname,'-dpdf');
    close(s)

    figure(f2)
    subplot(length(rois),2,(ri-1)*2+1); hold on
%     plot(ucon,mean(mat{1},2),'-','Color',clist1(end,:));
%     plot(ucon,mean(mat{2},2),'-','Color',clist2(end,:));
%     legend({'Cued Coh','Cued Con'});
    % get the pooled se...
    pooledse{1} = sum(matse{1}.*(Nm{1}-1),2) ./ sum(Nm{1}-1,2);
    pooledse{2} = sum(matse{2}.*(Nm{2}-1),2) ./ sum(Nm{1}-1,2);
    % pooled
    errorbar(ucon,mean(mat{1},2),1.96*pooledse{1},'o','Color',clist1(end,:));
    % regression
    b = [ones(length(ucon),1) ucon']\mean(mat{1},2);
    yCalc = b(1) + b(2) * ucon;
    plot(ucon,yCalc,'-','Color',clist1(end,:));
    %next
    errorbar(ucon,mean(mat{2},2),1.96*pooledse{2},'o','Color',clist2(end,:));
    % regression
    b = [ones(length(ucon),1) ucon']\mean(mat{2},2);
    yCalc = b(1) + b(2) * ucon;
    plot(ucon,yCalc,'-','Color',clist2(end,:));
    %next
    title(sprintf('%s Orange = Attend Con; Blue = Attend Coh',roi));
    xlabel('Contrast');
    ylabel('Response Amplitude');
    miv = min([mean(mat{1},2)-1.96.*pooledse{1};mean(mat{2},2)-1.96.*pooledse{2}])-.05;
    mav = max([mean(mat{1},2)+1.96.*pooledse{1};mean(mat{2},2)+1.96.*pooledse{2}])+.05;
    axis([min(ucon)-.05 max(ucon)+.05 miv mav]);
    
    % REPEAT FOR SP
    [s,f] = sp_figure(); hold on
    % pooled
    errorbar(ucon,mean(mat{1},2),1.96*pooledse{1},'o','Color',clist1(end,:));
    % regression
    b = [ones(length(ucon),1) ucon']\mean(mat{1},2);
    yCalc = b(1) + b(2) * ucon;
    plot(ucon,yCalc,'-','Color',clist1(end,:));
    roi_effects(ri,2,1,:) = yCalc;
    %next
    errorbar(ucon,mean(mat{2},2),1.96*pooledse{2},'o','Color',clist2(end,:));
    % regression
    b = [ones(length(ucon),1) ucon']\mean(mat{2},2);
    yCalc = b(1) + b(2) * ucon;
    plot(ucon,yCalc,'-','Color',clist2(end,:));
    roi_effects(ri,2,2,:) = yCalc;
    %next
    title(sprintf('%s Orange = Attend Con; Blue = Attend Coh',roi));
    xlabel('Contrast');
    ylabel('Response Amplitude');
    axis([min(ucon)-.05 max(ucon)+.05 miv mav]);
    fname = fullfile(dir1,'naka',sprintf('mean_contrast_%s',roi));
    print(fname,'-dpdf');
    close(s)

    subplot(length(rois),2,(ri-1)*2+2); hold on
%     plot(ucoh,mean(mat{1},1),'o','Color',clist1(end,:));
%     plot(ucoh,mean(mat{2},1),'o','Color',clist2(end,:));
%     legend({'Cued Coh','Cued Con'});
    % get the pooled se...
    pooledse{1} = sum(matse{1}.*(Nm{1}-1),1) ./ sum(Nm{1}-1,1);
    pooledse{2} = sum(matse{2}.*(Nm{2}-1),1) ./ sum(Nm{1}-1,1);
    % pooled
    errorbar(ucoh,mean(mat{1},1),1.96*pooledse{1},'o','Color',clist1(end,:));    
    % regression
    b = [ones(length(ucoh),1) ucoh']\mean(mat{1},1)';
    yCalc = b(1) + b(2) * ucoh;
    roi_effects(ri,1,1,:) = yCalc;
    plot(ucoh,yCalc,'-','Color',clist1(end,:));
    %next
    errorbar(ucoh,mean(mat{2},1),1.96*pooledse{2},'o','Color',clist2(end,:));
    % regression
    b = [ones(length(ucoh),1) ucoh']\mean(mat{2},1)';
    yCalc = b(1) + b(2) * ucoh;
    roi_effects(ri,1,2,:) = yCalc;
    plot(ucoh,yCalc,'-','Color',clist2(end,:));
    %next
    title(sprintf('%s Orange = Attend Con; Blue = Attend Coh',roi));
    xlabel('Coherence');
    ylabel('Response Amplitude');
    axis([min(ucoh)-.05 max(ucoh)+.05 miv mav]);
    
    % REPEAT FOR SP
    [s,f] = sp_figure(); hold on
   
    % pooled
    errorbar(ucoh,mean(mat{1},1),1.96*pooledse{1},'o','Color',clist1(end,:));    
    % regression
    b = [ones(length(ucoh),1) ucoh']\mean(mat{1},1)';
    yCalc = b(1) + b(2) * ucoh;
    plot(ucoh,yCalc,'-','Color',clist1(end,:));
    %next
    errorbar(ucoh,mean(mat{2},1),1.96*pooledse{2},'o','Color',clist2(end,:));
    % regression
    b = [ones(length(ucoh),1) ucoh']\mean(mat{2},1)';
    yCalc = b(1) + b(2) * ucoh;
    plot(ucoh,yCalc,'-','Color',clist2(end,:));
    %next
    title(sprintf('%s Orange = Attend Con; Blue = Attend Coh',roi));
    xlabel('Coherence');
    ylabel('Response Amplitude');
    axis([min(ucoh)-.05 max(ucoh)+.05 miv mav]);
    fname = fullfile(dir1,'naka',sprintf('mean_coherence_%s',roi));
    print(fname,'-dpdf');
    close(s)

end   
figure(f1);
    fname = fullfile(dir,'AmplitudeRange');
    print(fname,'-dpdf');
figure(f2);
fname = fullfile(dir,'AmplitudeMean');
print(fname,'-dpdf');
close(f1)
close(f2)
end

%% ROI plots
features = {'coherence','contrast'};
tasks = features;
figure('name',sprintf('ROI: %s',roi));
legvals{1} = {}; legvals{2} = {};
colors1 = brewermap(length(rois)*2,'Spectral');
for ri = 1:length(rois)
    [conVals, cohVals, cuedTask] = parseNames(neural.SCM.(name).(rois{ri}).stimNames); 
    ucoh = unique(cohVals);
    ucon = unique(conVals);
    vals = {ucoh,ucon};
    roi = rois{ri};
    for feature = 1:2
        subplot(1,2,feature), hold on
        for task = 1:2
%             title(sprintf('Feature: %s, Task: %s',features{feature},tasks{task}));
            plot(vals{feature},squeeze(roi_effects(ri,feature,task,:)),'Color',colors1((task-1)*length(rois)+ri,:));
            legvals{feature}{end+1} = sprintf('%s, attend %s',roi,tasks{task});
        end
        subplot(1,2,feature)
        xlabel(features{feature})
        legend(legvals{feature})
    end
end
fname = fullfile(dir,'ROI');
print(fname,'-dpdf');
close(gcf)

function [f,c] = sp_figure()

c = gcf;
f = figure;
whitebg(f);

function [conVal, cohVal, cuedTask] = parseNames(stimNames)
%% Note that parse names ignores the prefix entirely, so it's hemisphere independent
conVal = []; cohVal = []; cuedTask = [];

for i = 1:length(stimNames)
    name = stimNames{i};
    ands = strfind(name,' and');
    
    if strfind(name,'Con=')
        conVal(end+1) = str2num(name(strfind(name,'Con=')+4:ands(1))); % note contrast first
    end
    if strfind(name,'Coh=')
        cohVal(end+1) = str2num(name(strfind(name,'Coh=')+4:ands(2)-2));
    end
    % t can be 1, 2, (coherence, contrast) or 3, 4 (catch coherence, catch
    % contrast)
    cuedTask = [cuedTask str2num(name(strfind(name,'ask=')+4:end))];
end

