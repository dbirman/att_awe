function cc_eyePlots(sid)

allData = loadAllData(sid);

eye = allData.behav.eye;

%% Plot first 500 ms for all trials, overlapped

f500 = figure;
hold on

mx= [];
my = [];
for ei = 1:length(eye)
    ceye = eye{ei};
    for ti = 1:min(5,size(ceye.xPos,1))
        mx = [mx ceye.xPos(1:250)];
        my = [my ceye.yPos(1:250)];
        plot(ceye.xPos(1:250),ceye.yPos(1:250),'*');
    end
end
            mx = nanmean(mx); my = nanmean(my);

rectangle('Position',[mx-.5 my-.5 1 1],'Curvature',[1 1]);
rectangle('Position',[mx+3.5 my-5 7.5 10]);
rectangle('Position',[mx+-3.5-7.5 my-5 7.5 10]);