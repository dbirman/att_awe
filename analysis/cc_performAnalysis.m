function [neural, CRF] = cc_performAnalysis(neural,name)

neural.fit = struct;

rois = neural.shortROIs;

[conVals, cohVals, cuedTask] = parseNames(neural.SCM.(name).(roi).stimNames); 
%% Calculate Non-Lin Fits (for amplitude)

% fitter = 'fitType=nonlin';
% amper = 'amplitudeType=fit2';
% fits = {};
% for ri = 1:length(rois)
%     roi = rois{ri};       
%     
%     fits{end+1} = fitTimecourse(neural.tSeries.(name).(roi).tSeries,neural.SCM.(name).(roi).stimVol,.5,'concatInfo',neural.tSeries.(name).(roi).concatInfo,fitter,amper);
% end

%% Calculate Deconvolution by Condition

% fitter = 'fitType=deconv';
% amper = 'amplitudeType=fit1';
% deconv = {};
% for ri = 1:length(rois)
%     roi = rois{ri};       
%     
%     [conVals, cohVals, cuedTask] = parseNames(neural.SCM.(name).(roi).stimNames); 
%     deconv{end+1} = fitTimecourse(neural.tSeries.(name).(roi).tSeries,neural.SCM.(name).(roi).stimVol,.5,'concatInfo',neural.tSeries.(name).(roi).concatInfo,fitter,amper);
% end

%     d.deconv = getr2timecourse(timecourse,d.nhdr,d.hdrlenTR,d.scm,d.framePeriod,d.verbose);
    deconv = {};
for ri = 1:length(rois)
    roi = rois{ri};       
    
    d = constructD(neural.tSeries.(name).(roi).tSeries,neural.SCM.(name).(roi).stimVol,0.5,20,neural.tSeries.(name).(roi).concatInfo,'none','deconv',0);
%     deconv{end+1} = getr2timecourse(,,.5,'concatInfo',,fitter,amper);
end



%% Deconv Figures

% Contrast Response at each Coherence
ucoh = unique(cohVals);
ucon = unique(conVals);
clist1 = brewermap(4,'Oranges');
clist2 = brewermap(4,'Purples');
if ~length(ucoh)==length(ucon)
    disp('This code will break.. fix it up');keyboard
end
tasks = {'Cued Coherence','Cued Contrast'};
for ri = 1:length(rois)
    roi = rois{ri};
    dfs(1) = figure('name',sprintf('%s: %s',tasks{1},roi)); dfs(2) = figure('name',sprintf('%s: %s',tasks{2},roi));
%     fit = fits{ri}.fit;
    fit = deconv{ri};
    
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
                plot(ehdr,'Color',clist1(ci,:));
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
                plot(ehdr,'Color',clist2(ci,:));
            end
            title(sprintf('Contrast: %0.02f',ccon));
        end
    end
end

%% Amplitude Figures
f1 = figure;
f2 = figure;
for ri = 1:length(rois)
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
    legend(legVals);
    for i = 1:4
        errorbar(ucon,mat{1}(:,i),matse{1}(:,i),'*','Color',clist1(i,:));
        errorbar(ucon,mat{2}(:,i),matse{2}(:,i),'*','Color',clist2(i,:));
    end
    title(sprintf('%s: Orange = Cued Coh, Green = Cued Con.',roi));
    xlabel('Contrast');
    ylabel('Response Amplitude');
    axis([0 1 0 ceil(max([mat{1}(:,i)',mat{2}(:,i)']))]);

    % X axis = coherence
    subplot(length(rois),2,(ri-1)*2+2); hold on
    for i = 1:4
        plot(ucoh,mat{1}(i,:),'-','Color',clist1(i,:));
        errorbar(ucoh,mat{1}(i,:),matse{1}(i,:),'*','Color',clist1(i,:));
        plot(ucoh,mat{2}(i,:),'-','Color',clist2(i,:));
        errorbar(ucoh,mat{2}(i,:),matse{2}(i,:),'*','Color',clist2(i,:));
    end
    title(sprintf('%s: Orange = Cued Coh, Green = Cued Con.',roi));
    xlabel('Coherence');
    ylabel('Response Amplitude');
    axis([0 1 0 ceil(max([mat{1}(i,:),mat{2}(i,:)]))]);


    figure(f2)
    subplot(length(rois),2,(ri-1)*2+1); hold on
    plot(ucon,mean(mat{1},2),'-','Color',clist1(end,:));
    plot(ucon,mean(mat{2},2),'-','Color',clist2(end,:));
    legend({'Cued Coh','Cued Con'});
    % get the pooled se...
    pooledse{1} = sum(matse{1}.*(Nm{1}-1),2) ./ sum(Nm{1}-1,2);
    pooledse{2} = sum(matse{2}.*(Nm{2}-1),2) ./ sum(Nm{1}-1,2);
    % pooled
    errorbar(ucon,mean(mat{1},2),1.96*pooledse{1},'*','Color',clist1(end,:));
    errorbar(ucon,mean(mat{2},2),1.96*pooledse{2},'*','Color',clist2(end,:));
    title(sprintf('%s: Orange = Cued Coh, Green = Cued Con.',roi));
    xlabel('Contrast');
    ylabel('Response Amplitude');
    axis([0 1 0 ceil(max(max([mean(mat{1},2),mean(mat{2},2)])))]);

    subplot(length(rois),2,(ri-1)*2+2); hold on
    plot(ucoh,mean(mat{1},1),'-','Color',clist1(end,:));
    plot(ucoh,mean(mat{2},1),'-','Color',clist2(end,:));
    legend({'Cued Coh','Cued Con'});
    % get the pooled se...
    pooledse{1} = sum(matse{1}.*(Nm{1}-1),1) ./ sum(Nm{1}-1,1);
    pooledse{2} = sum(matse{2}.*(Nm{2}-1),1) ./ sum(Nm{1}-1,1);
    % pooled
    errorbar(ucoh,mean(mat{1},1),1.96*pooledse{1},'*','Color',clist1(end,:));
    errorbar(ucoh,mean(mat{2},1),1.96*pooledse{2},'*','Color',clist2(end,:));
    title(sprintf('%s: Orange = Cued Coh, Green = Cued Con.',roi));
    xlabel('Coherence');
    ylabel('Response Amplitude');
    axis([0 1 0 ceil(max(max([mean(mat{1},1),mean(mat{2},1)])))]);
end   
    
%% TODO: Add print functions to print out figures
stop = 1;

%% CRF
for ri = 1:length(rois)
    roi = rois{ri};
    
    fdata = allData.fit.(roi);
    
    CRF.(sprintf('%s_%s',roi,'contrast')).cued.c = ucon;
    CRF.(sprintf('%s_%s',roi,'coherence')).cued.c = ucoh;
    
    CRF.(sprintf('%s_%s',roi,'contrast')).cued.r = mean(fdata.mat{2},2); % contrast pedestals, cued to contrast
    CRF.(sprintf('%s_%s',roi,'coherence')).cued.r = mean(fdata.mat{1},1); % coherence pedestals, cued to coherence
    
    
    CRF.(sprintf('%s_%s',roi,'contrast')).uncued.c = ucon;
    CRF.(sprintf('%s_%s',roi,'coherence')).uncued.c = ucoh;
    
    CRF.(sprintf('%s_%s',roi,'contrast')).uncued.r = mean(fdata.mat{1},2); % contrast pedestals, cued to contrast
    CRF.(sprintf('%s_%s',roi,'coherence')).uncued.r = mean(fdata.mat{2},1); % coherence pedestals, cued to coherence
end


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

























%%%%%%%%%%%%%%%%%%%%
%%   constructD   %%
%%%%%%%%%%%%%%%%%%%%
function d = constructD(timecourse,stimvol,framePeriod,hdrlen,concatInfo,option,fitType,verbose)

% make sure timecourse looks ok.
if (size(timecourse,2) == 1) && (size(timecourse,1) > 1)
  mrWarnDlg(sprintf('(fitTimecourse) Timecourse found to be %ix1 instead of 1x%i. Taking transpose',size(timecourse,1),size(timecourse,1)));
  timecourse = timecourse';
end

% first get stimulus convolution matrix for each condition separately
d.dim(4) = size(timecourse,2);
% note that this hdrlen is different from
% below, becuase it is in # of TR's not
% in seconds.
hdrlenTR = floor(hdrlen/framePeriod)+1;
d.hdrlen = hdrlenTR;
d.concatInfo = concatInfo;
d.nFrames = d.dim(4);

if verbose
  disppercent(-inf,'(fitTimecourse) Constructing stimulus convolution matrices for each condition');
end
for stimNum = 1:length(stimvol)
  d.stimvol{1} = stimvol{stimNum};
  d = makescm(d,[],0);
  eachSCM{stimNum} = d.scm;
  if verbose,disppercent(stimNum/length(stimvol));end
end
if verbose,disppercent(inf);end

% make an scm for the whole design
if verbose
  disppercent(-inf,'(fitTimecourse) Constructing stimulus convolution matrix for whole design');
end
d.stimvol = stimvol;
d = makescm(d,[],1);
allSCM = d.scm;
if verbose,disppercent(inf);end

% make an scm for each stimulus 
% now make a d structure to pass around
clear d;
d.dim = size(timecourse);
d.stimvol = stimvol;
d.framePeriod = framePeriod;
d.concatInfo = concatInfo;
d.eachSCM = eachSCM;
d.scm = allSCM;
d.hdrlen = hdrlen;
d.hdrlenTR = hdrlenTR;
d.nhdr = length(stimvol);
d.timecourse = timecourse;
d.applyFiltering = 1;
d.zeroMean = 0;
d.deconvModel = 1;
d.nFrames = max(d.dim);
d.verbose = verbose;