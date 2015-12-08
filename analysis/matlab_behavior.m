%% Setup


% subjects = {'s021','s025','s300','s302','s304','s305','s309','s310'};
subjects = {'s300','s304','s305','s315'};

%% Copy Files

for si = 1:length(subjects)
    sid = subjects{si};
    
    cc_copyFiles(sid);
end

%% Load and Export

allData = loadAllData(sid);

maf = {};
caf = {};
nocaf = {};
for si = 1:length(subjects)
    sid = subjects{si};
    
    mglSetSID(sid);
    
    allData = exportCohCon(allData);
    saveAllData(sid,allData);
end


%% Eye Plots

for si = 1:length(subjects)
    sid = subjects{si};
    
    cc_eyePlots(sid);
end

%% R-Choice Plots

cc_rightChoicePlots(true);   

%% State-Trace Plots

for si = 1:length(subjects)
    sid = subjects{si};
        
    [f_stp.(sid).f, f_stp.(sid).data] = cc_stateTracePlot(sid);
end

%% Generate fullData
% fullData is the entire dataset in long form, it loads all of the subjects
% individually and copies their entire dataset (every trial) into a single
% matrix. The formatting is as follows:
% Run # | Trial # | lCon | rCon | lCoh | rCoh | dCon | dCoh | side | response | task 

fullData = loadFullData();
reset = 1;

if reset==1 || length(fields(fullData))==0 || isempty(fullData.data) || ~length(subjects)==length(fullData.subjects)
    fullData.data = []; fullData.header = {}; fullData.subjects = {};
    fullData.subjects = subjects;
    fData = {};
    for si = 1:length(subjects)
        sid = subjects{si};

        [fData{si}, header] = cc_exportFull(sid);
        fData{si}(:,end+1) = mrStr2num(strrep(sid,'s',''));
    end
    header{end+1} = 'sid';

    data = cc_catLong(fData);

    % add 'side'
    header{end+1} = 'side';

    side = [];
    for i = 1:size(data)
        if data(i,3)==1 % coherence
            side(i) = data(i,10);
        else
            side(i) = data(i,9);
        end
    end
    data(:,end+1) = side;

    fullData.data = data;
    fullData.header = header;
end

saveFullData(fullData);

%% lightning talk plot

clist2 = brewermap(3,'Oranges');
clist1 = brewermap(3,'Blues');

mafs = cell(length(subjects),1); cafs = cell(length(subjects),1);
for si = 1:length(subjects)
    sid = subjects{si};
    allData = loadAllData(sid);
    
    figure
    subplot(1,2,1), hold on
    % attending coherence
    maf1 = allData.behav.mfits{1}(:,1);
    caf1 = allData.behav.cfits{1}(:,1);
    bar([1 2],[mean(maf1) mean(caf1)]);
    errorbar([1 2],[mean(maf1) mean(caf1)],[std(maf1) std(caf1)],'*');
    axis([0.5 2.5 0 max([mean(maf1) mean(caf1)]+[std(maf1) std(caf1)])]);
    subplot(1,2,2), hold on
    % attending contrast
    maf2 = allData.behav.mfits{2}(:,1);
    caf2 = allData.behav.cfits{2}(:,1);
    b = bar([1 2],[mean(maf2) mean(caf2)]);
    errorbar([1 2],[mean(maf2) mean(caf2)],[std(maf2) std(caf2)],'*');
    axis([0.5 2.5 0 max([mean(maf2) mean(caf2)]+[std(maf2) std(caf2)])]);
    dir = fullfile('~/proj/att_awe/analysis/figures',sid);
    if ~isdir(dir), mkdir(dir); end
    fname = fullfile(dir,'threshold_bar.pdf');
    print(fname,'-dpdf');
    
    mafs{si,1} = maf1;
    cafs{si,1} = caf1;
    mafs{si,2} = maf2;
    cafs{si,2} = caf2;
end

%% average plot

marray = cellfun(@mean, mafs);
carray = cellfun(@mean, cafs);

mu = mean(marray,1);
ms = std(marray,1);
cu = mean(carray,1);
cs = std(carray,1);

topy = [0.55 0.21];

figure
for pos = 1:2
    subplot(1,2,pos), hold on
    bar([1 2],[mu(pos) cu(pos)]);
    errorbar([1 2],[mu(pos) cu(pos)],[ms(pos) cs(pos)],'o');
    axis([0.5 2.5 0 topy(pos)]);
    fname=fullfile('~/proj/att_awe/analysis/figures','avg_bar.pdf');
    print(fname,'-dpdf');
    [h,p] = ttest2(marray(:,pos),carray(:,pos)); p
end