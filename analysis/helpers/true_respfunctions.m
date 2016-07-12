function true_respfunctions(subjs,mode,rois)
%UNTITLED2 Generate a plot of the 'true' response functions (including
%offsets in both attention conditions


if strfind(mode,'subj')
    for si = 1:length(subjs)
        plotsubj(subjs{si},rois);
    end
end

if strfind(mode,'all')
    plotall(subjs,rois);
end

function plotsubj(subj,rois)
%%
load(fullfile(sprintf('C:/Users/Dan/proj/COHCON_DATA/%s_fitatt.mat',subj)));

h = figure;
xlabel('Conrast/Coherence (%)');
ylabel('Response (% Signal Change / sec)');
% plot contrast neutral
cmap = brewermap(3,'PuOr');
x = 0:.01:1;
% stage info
sub = 1:length(rois);
subs = length(rois);
col = [1 3];
lines = {'-','--','-.'};
models = {'con','coh'};
modelFunc = {'conModel','cohModel'};
attnames = {'neutral','conatt','cohatt'};
attvals = [0 2 1];

for ri = 1:length(rois)
    params.(rois{ri}) = mergeroiparams(fitatt.roiparams,find(cell2mat(cellfun(@(x) ~isempty(strfind(x,rois{ri})),fitatt.ROIs,'UniformOutput',false))));
end
for i = 1:length(rois)
    subplot(subs,1,sub(i)), hold on
    for k = 1:2
        for j = 1:3
            % get data
            eval(sprintf('temp = %s(x,params.%s,%i,0);',modelFunc{k},rois{i},attvals(j)));
            eval(sprintf('plot(x,temp,''%s'',''Color'',cmap(%i,:));',lines{j},col(k)));
        end
    end
    xlabel('Contrast/Coherence (%)');
    ylabel('Response (% Signal Change / sec)');
    title(rois{i});
    drawPublishAxis
end
%run

% fname = fullfile('C:/Users/Dan/proj/COHCON_DATA',sprintf('%s_fullroiresp.pdf',subj));
% set(h,'Units','Inches');
% pos = get(h,'Position');
% set(h,'InvertHardCopy','off');
% set(gcf,'Color',[1 1 1]);
% set(gca,'Color',[1 1 1]);
% set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(fname,'-dpdf');

function plotall(subjs,rois)
%%
%we'll do the simplest thing, just average the roiparams for all the
%subjects together, and then plot as above 
bigroiparams = [];
locs = cell(1,20);
for i = 1:length(subjs)
    subj = subjs{i};
    load(fullfile(sprintf('C:/Users/Dan/proj/COHCON_DATA/%s_fitatt.mat',subj)));
    %fitatt.roiparams containts data
    bigroiparams = [bigroiparams fitatt.roiparams];
    for j = 1:20
        locs{j} = [locs{j} j+(i-1)*20];
    end
end
ROIs = fitatt.ROIs;
fitatt = struct;
fitatt.ROIs = ROIs;
for i = 1:20
    fitatt.roiparams{i} = mergeroiparams(bigroiparams,locs{i});
end
save(fullfile('C:/Users/Dan/proj/COHCON_DATA/avg_fitatt.mat'),'fitatt');
plotsubj('avg',rois);