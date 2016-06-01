function thresholdTimePlot( stimulus )
%thTimePlot Make plots of threshold over time for cohcon

%%
stair = stimulus.staircases;

% main
figure, hold on

main = stair.main;
clist = brewermap(3,'PuOr');
convals = [];
for i = 1:length(main{2}) 
    try
    out = doStaircase('threshold',main{2}(i));
    catch
        out.threshold = NaN;
    end
    convals = [convals out.threshold];
end
cohvals = [];
for i = 1:length(main{1})
    try
        out = doStaircase('threshold',main{1}(i));
    catch
        out.threshold = NaN;
    end
    cohvals = [cohvals out.threshold];
end

%%
clf
hold on
plot(convals,'-o','MarkerFaceColor',clist(1,:),'MarkerEdgeColor',[1 1 1]);
plot(cohvals,'-o','MarkerFaceColor',clist(3,:),'MarkerEdgeColor',[1 1 1]);
