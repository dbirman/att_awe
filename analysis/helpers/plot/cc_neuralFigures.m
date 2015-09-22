function cc_neuralFigures(neural,name,sid,df, nf)

dir = fullfile('~/proj/att_awe/analysis/figures',sid);
if ~isdir(dir), mkdir(dir); end

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
        figure(dfs(task));
        title(sprintf('Task %s',tasks{task}));
        for ui = 1:length(ucoh)
            ccoh = ucoh(ui);
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
                plot(ehdr,'Color',clist1(ci,:));
%                 errorbar(ehdr,se,'Color',clist1(ci,:));
            end
            title(sprintf('Coherence: %0.02f',ccoh));
        end
        for ui = 1:length(ucon)
            ccon = ucon(ui);
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
                plot(ehdr,'Color',clist2(ci,:));
%                 errorbar(ehdr,se,'Color',clist2(ci,:));
            end
            title(sprintf('Contrast: %0.02f',ccon));
        end
        
        fname = fullfile(dir,sprintf('Deconv_%s',roi));
        print(fname,'-dpdf');

    end
end
end

%% Amplitude Figures
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
    % X axis = contrast
    subplot(length(rois),2,(ri-1)*2+1); hold on
    legVals = {};
    for i = 1:4
        plot(ucon,mat{1}(:,i),'-','Color',clist1(i,:));
        legVals{end+1} = sprintf('Cued Coh, Coh = %0.2f',ucoh(i));
    end
    for i = 1:4
        plot(ucon,mat{2}(:,i),'-','Color',clist2(i,:));
        legVals{end+1} = sprintf('Cued Con, Coh = %0.2f',ucoh(i));
    end
%     legend(legVals);
    for i = 1:4
        errorbar(ucon,mat{1}(:,i),matse{1}(:,i),'*','Color',clist1(i,:));
        errorbar(ucon,mat{2}(:,i),matse{2}(:,i),'*','Color',clist2(i,:));
    end
    title(sprintf('%s Orange = Attend Con; Blue = Attend Coh',roi));
    xlabel('Contrast');
    ylabel('Response Amplitude');
    axis([min(ucon)-.05 max(ucon)+.05 round(10*min([mat{1}(i,:),mat{2}(i,:)]))/10-.1 round(10*max([mat{1}(:,i)',mat{2}(:,i)']))/10+.1]);
%     set(ax,'YTick',[-1 -0.5 0 0.5 1])
    
    % X axis = coherence
    subplot(length(rois),2,(ri-1)*2+2); hold on
    for i = 1:4
        plot(ucoh,mat{1}(i,:),'-','Color',clist1(i,:));
        errorbar(ucoh,mat{1}(i,:),matse{1}(i,:),'*','Color',clist1(i,:));
        plot(ucoh,mat{2}(i,:),'-','Color',clist2(i,:));
        errorbar(ucoh,mat{2}(i,:),matse{2}(i,:),'*','Color',clist2(i,:));
    end
    title(sprintf('%s Orange = Attend Con; Blue = Attend Coh',roi));
    xlabel('Coherence');
    ylabel('Response Amplitude');
    axis([min(ucoh)-.05 max(ucoh)+.05 round(10*min([mat{1}(i,:),mat{2}(i,:)]))/10-.1 round(10*max([mat{1}(i,:),mat{2}(i,:)]))/10+.1]);


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
    axis([min(ucon)-.05 max(ucon)+.05 round(10*min(min([mean(mat{1},2),mean(mat{2},2)])))/10-.1 round(10*max(max([mean(mat{1},2),mean(mat{2},2)])))/10+.1]);

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
    axis([min(ucoh)-.05 max(ucoh)+.05 round(10*min(min([mean(mat{1},1),mean(mat{2},1)])))/10-.1 round(10*max(max([mean(mat{1},1),mean(mat{2},1)])))/10+.1]);

end   
figure(f1);
    fname = fullfile(dir,'AmplitudeRange');
    print(fname,'-dpdf');
figure(f2);
fname = fullfile(dir,'AmplitudeMean');
print(fname,'-dpdf');
    
%% CRF
% for ri = 1:length(rois)
%     roi = rois{ri};
%     
%     fdata = fits{ri};
%     
%     CRF.(sprintf('%s_%s',roi,'contrast')).cued.c = ucon;
%     CRF.(sprintf('%s_%s',roi,'coherence')).cued.c = ucoh;
%     
%     CRF.(sprintf('%s_%s',roi,'contrast')).cued.r = mean(fdata.mat{2},2); % contrast pedestals, cued to contrast
%     CRF.(sprintf('%s_%s',roi,'coherence')).cued.r = mean(fdata.mat{1},1); % coherence pedestals, cued to coherence
%     
%     
%     CRF.(sprintf('%s_%s',roi,'contrast')).uncued.c = ucon;
%     CRF.(sprintf('%s_%s',roi,'coherence')).uncued.c = ucoh;
%     
%     CRF.(sprintf('%s_%s',roi,'contrast')).uncued.r = mean(fdata.mat{1},2); % contrast pedestals, cued to contrast
%     CRF.(sprintf('%s_%s',roi,'coherence')).uncued.r = mean(fdata.mat{2},1); % coherence pedestals, cued to coherence
% end



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

