%% Load the png with ROIs

png = imread('~/proj/att_awe/talks/lightning_2018/base_withrois.png');

%% find the v1 area
for x=1:361
    for y=1:361
        v1(x,y) = png(x,y,1)==153 && png(x,y,2)==142 && png(x,y,3)==195;
        mt(x,y) = png(x,y,1)==247 && png(x,y,2)==169 && png(x,y,3)==0;
    end
end



%% load the data
load /Users/dan/proj/att_awe/analysis/osf_predict_func/data.mat

%% grab the v1 and mt 100% response functions
con = data{1}.adata{1}.cc.con;
coh = data{1}.adata{1}.cc.coh;

conidx = 16;
cohidx = 5;
maxidx = 20;

rois = [1 8];
idxs = [conidx cohidx maxidx];

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

%% 

figure;
plot(squeeze(max_(:,3,:))');

%% For each index, generate a gif
map = brewermap(100,'Oranges');
range = 0:3/99:3;

figure;
for ii= 1:length(idxs)
    gif(fullfile('~/proj/att_awe/talks/lightning_2018/gifs/',sprintf('%i.gif',idxs(ii))),'DelayTime',0.500);
    
    % get the data
    v1dat = squeeze(max_(1,ii,:));
    mtdat = squeeze(max_(2,ii,:));
    
    r = png(:,:,1);
    g = png(:,:,2);
    b = png(:,:,3);
    
    for ti = 1:length(v1dat)
        rv1 = v1dat(ti);
        rmt = mtdat(ti);
        % separate the channels
        r(v1) = interp1(range,map(:,1),rv1);
        g(v1) = interp1(range,map(:,2),rv1);
        b(v1) = interp1(range,map(:,3),rv1);
        
        r(mt) = interp1(range,map(:,1),rmt);
        g(mt) = interp1(range,map(:,2),rmt);
        b(mt) = interp1(range,map(:,3),rmt);
        
        tpng(:,:,1) = r;
        tpng(:,:,2) = g;
        tpng(:,:,3) = b;
        
        imagesc(tpng);
        pause(.01);
        
        gif;
    end
end