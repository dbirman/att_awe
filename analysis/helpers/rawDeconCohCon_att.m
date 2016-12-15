function rawDeconCohCon_att( roiname , subj)
%DECONCOHCON Deconvolve the cohxcon experiment

%% define folders
datafiles = dir(fullfile(datafolder,sprintf('%s*_data.mat',subj)));

%% load data
datas = {};
for fi = 1:length(datafiles)
    load(fullfile(datafolder,datafiles(fi).name));
    % we only want to keep non-fixation task, so check >0
    if all(data.design(:,9)>0)
        datas{end+1}=data;
    end
end

data = datas;

if length(data)<2
    warning(sprintf('Only %i files found for %s',length(data),subj));
    return
end

keep = ones(size(data));
for di = 1:length(data)
    if isempty(data{di}.tSeries)
        keep(di) = 0;
    end
end
data = data(logical(keep));
%% Concatenate Sessions if Necessary

tSeriesname = 'rtSeries25';
if iscell(data) && length(data)>1
    disp('(roimodel) Concatenating what appear to be different sessions...');
    % to run the model we really only need:
    % data.tSeries
    % data.design
    % data.runtrans
    % but we need these to be concatenated across the different sessions.
    % We have to be a little careful here that we make sure we add the
    % lengths of each run correctly or we'll screw ourselves over.
    data_old = data;
    data = struct;
    data.ROIs = data_old{1}.ROIs;
    % (we can do this the slow way, perf doesn't really matter
    data.tSeries = cell(1,length(data_old{1}.ROIs));
    length_sofar = 0;
    for ri = 1:length(data.tSeries), data.tSeries{ri} = []; end
    data.design = [];
    data.runtrans = [];
    for si = 1:length(data_old)
        for ri = 1:length(data_old{si}.tSeries)
            data.tSeries{ri} = [data.tSeries{ri} data_old{si}.(tSeriesname){ri}];
        end
        % tweak the SV by adding 
        cdes = data_old{si}.design;
        cdes(:,1) = cdes(:,1) + length_sofar;
        data.design = [data.design ; cdes];
        data.runtrans = [data.runtrans ; data_old{si}.runtrans+length_sofar];
        length_sofar = length_sofar + length(data_old{si}.tSeries{1});
    end
else
    data = data{1};
end

%%
if isempty(strfind(roiname,'l')) && isempty(strfind(roiname,'r'))
    warning('rawDecon has to be called with a single ROI because the responses are l/r specific!');
    keyboard
end

if isempty(strfind(roiname,'l'))
    conidx = 3;
    cohidx = 6;
elseif isempty(strfind(roiname,'r'))
    conidx = 4;
    cohidx = 7;
end

roinum = find(cellfun(@(x) strcmp(x,roiname),data.ROIs),1);
% DON'T CONCATENATE ROIS!!!

%% Rebuild the design
%    1        2    3     4     5      6   7    8      9
% stimvol basecon lcon rcon basecoh lcoh rcoh timing task

design = data.design;
%% Reduce space of cohcon values (rounding to nearest true shown value)

con_idxs = [3 4];
coh_idxs = [6 7];

conopts = [0.325 0.4 0.55 0.85];
cohopts = [0.15 0.3 0.45 0.6];

for i = 1:length(con_idxs)
    for j = 1:size(design,1)
        con = design(j,con_idxs(i));
        low = find(con>=conopts,1,'last');
        if low==length(conopts)
            ncon = conopts(end);
        elseif abs(con-conopts(low)) > abs(con-conopts(low+1))
            ncon = conopts(low+1);
        else
            ncon = conopts(low);
        end
        design(j,con_idxs(i)) = ncon;
        
        coh = design(j,coh_idxs(i));
        low = find(coh>=cohopts,1,'last');
        if low==length(cohopts)
            ncoh = cohopts(end);
        elseif abs(coh-cohopts(low)) > abs(coh-cohopts(low+1))
            ncoh = cohopts(low+1);
        else 
            ncoh = cohopts(low);
        end
        design(j,coh_idxs(i)) = ncoh;
    end
end

%% CohxCon
contrast = design(:,conidx); ucon = unique(contrast);
coherence = design(:,cohidx); ucoh = unique(coherence);
task = design(:,9); utask = unique(task);
sv = design(:,1);

cohxcon_sv = {};
% stimNames = {};
conidx = [];
cohidx = [];
taskidx = []; % split by task
for coni = 1:length(ucon)
    for cohi = 1:length(ucoh)
        for ti = 1:length(utask)
            cohxcon_sv{end+1} = sv(logical((task==utask(ti)).*(contrast==ucon(coni)).*(coherence==ucoh(cohi))));
%             stimNames{end+1} = sprintf('Contrast=%0.2f Coherence=%0.2f',ucon(coni),ucoh(cohi));
            conidx(end+1) = ucon(coni);
            cohidx(end+1) = ucoh(cohi);
            taskidx(end+1) = utask(ti);
        end
    end
end

concatInfo.runTransition = data.runtrans;
curd = constructD(data.tSeries{roinum},cohxcon_sv,0.5,40,concatInfo,'none','deconv',0);
decon = getr2timecourse(curd.timecourse,curd.nhdr,curd.hdrlenTR,curd.scm,curd.framePeriod,curd.verbose);
decon = rmfield(decon,'scm');
decon = rmfield(decon,'covar');

%% save cohxcon data
fname = fullfile(datafolder,sprintf('%s_decon_att.mat',subj));
if exist(fname,'file')==2, load(fname); end
decondata.(roiname).cc.conidxs = conidx;
decondata.(roiname).cc.cohidxs = cohidx;
decondata.(roiname).cc.task = taskidx;
decondata.(roiname).cc.resp = decon.ehdr;
save(fname,'decondata');
return
%%

clist = brewermap(10,'PuOr');
f= figure;
% plot #1: contrast against 0% coherence (no change)
lconidx = find(cohidx==0.15);
convalues = [0.325 0.4 0.55 0.85];
colmap = [4 3 2 1];
subplot(2,1,1), hold on
for i = 1:4
    ci = colmap(find(convalues==conidx(lconidx(i)),1));
    plot(decon.time,decon.ehdr(lconidx(i),:),'o','MarkerSize',10,'MarkerFaceColor',clist(ci,:),'MarkerEdgeColor',[1 1 1]);
    errbar(decon.time,decon.ehdr(lconidx(i),:),decon.ehdrste(lconidx(i),:),'Color',clist(ci,:));
end
for i = 1:4
    ci = colmap(find(convalues==conidx(lconidx(i)),1));
    plot(model.time,model.ehdr(lconidx(i),:),'Color',clist(ci,:));
end
% legend({'Contrast=+0%','Contrast=+25%','Contrast=+50%','Contrast=+75%'});
% plot #2: coherence against 25% contrast (no change)
lcohidx = find(conidx==convalues(1));
a = axis;
axis([0 14 a(3) a(4)]);
drawPublishAxis
subplot(2,1,2), hold on
for i=1:5
    plot(decon.time,decon.ehdr(lcohidx(i),:),'o','MarkerSize',10,'MarkerFaceColor',clist(i+5,:),'MarkerEdgeColor',[1 1 1]);
    errbar(decon.time,decon.ehdr(lcohidx(i),:),decon.ehdrste(lcohidx(i),:),'Color',clist(i+5,:));
end
for i = 1:5
    plot(model.time,model.ehdr(lcohidx(i),:),'Color',clist(i+5,:));
end
% legend({'Coherence+0%','Coherence+25%','Coherence+50%','Coherence+75%','Coherence+100%'});

a = axis;
axis([0 14 a(3) a(4)]);
drawPublishAxis
%% Print Figure #2

fname = fullfile(datafolder,sprintf('%s_%s_cohconplot.pdf',subj,roiname));
savepdf(f,fname);

%% Print Figure #1
title(sprintf('%s: %s',subj,roiname));
savepdf(f,fullfile(datafolder,sprintf('%s_%s_timeplot.pdf',subj,roiname)));
