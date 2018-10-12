%% Load the png with ROIs

png = imread('~/proj/att_awe/talks/lightning_2018/base_withrois.png');

%% find the v1 area
for x=1:361
    for y=1:361
        mt(x,y) = png(x,y,1)==153 && png(x,y,2)==142 && png(x,y,3)==195;
        v1(x,y) = png(x,y,1)==247 && png(x,y,2)==169 && png(x,y,3)==0;
    end
end

flatpng = png(:);
v1_ = repmat(v1,1,1,3);
flatv1 = v1_(:);
mt_ = repmat(mt,1,1,3);
flatmt = mt_(:);


%% load the data
load /Users/dan/proj/att_awe/analysis/osf_predict_func/data.mat

%% grab the v1 and mt 100% response functions
con = data{1}.adata{1}.cc.con;
coh = data{1}.adata{1}.cc.coh;

conidx = 16;
cohidx = 5;
maxidx = 20;

rois = [1 8];
idxs = [1:20];

clear max_resp

for si = 1:11
    for ri = 1:length(rois)
        roi = rois(ri);
        for ii = 1:length(idxs)
            idx = idxs(ii);
            max_resp(si,ri,ii,:) = squeeze(data{si}.adata{1}.cc.resp(roi,idx,:));
            max_resp(si,ri,ii,:) = squeeze(data{si}.adata{1}.cc.resp(roi,idx,:));
        end
    end
end

max_ = squeeze(mean(max_resp));

%% For each index, generate a gif
map = 255*brewermap(100,'Oranges');
time = 0.25:0.5:40.5;
rtime = 0.1:0.1:40.5;
range = fliplr(-1:4/99:3)';

h = figure; hold on
imshow(png);
for ii= 1:length(idxs)
    % get the data
    v1dat = squeeze(max_(1,ii,:));
    mtdat = squeeze(max_(2,ii,:));
    
    
    for ti = 1:250
        tpng = png;
        r = png(:,:,1);
        g = png(:,:,2);
        b = png(:,:,3);
        rv1 = interp1(time,v1dat,rtime(ti));
        if isnan(rv1)
            rv1 = v1dat(1);
        end
        rmt = interp1(time,mtdat,rtime(ti));
        if isnan(rmt)
            rmt = mtdat(1);
        end
        % separate the channels
        r(v1) = interp1(range,map(:,1),rv1);
        g(v1) = interp1(range,map(:,2),rv1);
        b(v1) = interp1(range,map(:,3),rv1);
        
        r(mt) = interp1(range,map(:,1),rmt);
        g(mt) = interp1(range,map(:,2),rmt);
        b(mt) = interp1(range,map(:,3),rmt);
        
        tpng(:,:,1) = reshape(r,size(v1));
        tpng(:,:,2) = reshape(g,size(v1));
        tpng(:,:,3) = reshape(b,size(v1));
        
        imshow(tpng);
        
        if ti==1
            gif(fullfile('~/proj/att_awe/talks/lightning_2018/gifs/',sprintf('%i.gif',idxs(ii))),'DelayTime',0.100,'LoopCount',1,'nodither');
        else
            gif('frame',gca);
        end
    end
end

%% Create a colorbar png
figure; hold on
map_ = reshape(map,100,1,3);
map_ = repmat(map_,1,20,1);
imshow(map_/255);
% for i = 1:100
%     
% end


%% Other plots
out = cohcon_predict({'V1','MT'},0:.01:1,0:.01:1,1);

h = figure;
plot(0:.01:1,out.contrastResponse(:,1),'Color',cmap(2,:),'LineWidth',2)
axis([0 1 0 2]);
set(gca,'XTick',[0 0.5 1],'XTickLabel',[0 50 100],'YTick',[0 1 2]);
xlabel('Contrast (%)');
ylabel('Signal change (%)');
drawPublishAxis('figSize=[6,4]');
savepdf(h,'~/proj/att_awe/talks/lightning_2018/con.pdf');

h = figure;
plot(0:.01:1,out.coherenceResponse(:,2),'Color',cmap(6,:),'LineWidth',2)
axis([0 1 0 2]);
set(gca,'XTick',[0 0.5 1],'XTickLabel',[0 50 100],'YTick',[0 1 2]);
xlabel('Coherence (%)');
ylabel('Signal change (%)');
drawPublishAxis('figSize=[6,4]');
savepdf(h,'~/proj/att_awe/talks/lightning_2018/coh.pdf');
