
%%%%%%%%%%%%%%%%%%%%%%%
%    dispInfo    %
%%%%%%%%%%%%%%%%%%%%%%%
function h = ccDispInfo(stimulus,subj)

%% No-Catch Performance
nocatch = zeros(2,4);
nocatchs = zeros(2,4);

for task = 1:size(nocatch,1)
    for ped = 1:size(nocatch,2)
        try
        out = doStaircase('threshold',stimulus.staircases.nocatch{task,ped},'type','weibull','dispFig=0');

%             out = doStaircase('threshold',stimulus.staircases.nocatch{task,ped});
        nocatch(task,ped) = out.threshold;
        nocatchs(task,ped) = out.thresholdSTE;
        catch % probably a missing staircase, whatever
        end
    end
end
if any(nocatch(:)>10)
    warning('removing stairs >10');
    nocatch(nocatch>10) = 0;
end
%% Main + Catch
main = zeros(2,1);
mains = zeros(2,1);
catch_ = zeros(2,1);
catch_s = zeros(2,1);

for task = 1:2
    try
        out = doStaircase('threshold',stimulus.staircases.main{task},'type','weibull','dispFig=0');
%             out = doStaircase('threshold',stimulus.staircases.main{task});
        main(task) = out.threshold;
        mains(task) = out.thresholdSTE;
    catch % missing a staircase
    end
    try 
        out = doStaircase('threshold',stimulus.staircases.catch{task},'type','weibull','dispFig=0');
%             out = doStaircase('threshold',stimulus.staircases.catch{task});
        catch_(task) = out.threshold;
        catch_s(task) = out.thresholdSTE;
    catch
    end
end

%% Save data
fname = fullfile(sprintf('C:/Users/Dan/Documents/Box Sync/COHCON_DATA/%s_data.mat',subj));
if exist(fname)==2, load(fname); end
data.control = nocatch;
data.attend = main;
data.unattend = catch_;
save(fname,'data');
    %% The Plot!
map = brewermap(6,'PuOr');

h = figure;
hold on
h2 = plot([0.15 0.3 0.45 0.6],nocatch(1,:),'-','Color',map(1,:));
%     errbar([0.15 0.3 0.45 0.6],nocatch(1,:),nocatchs(1,:),'Color',map(1,:));
h6 = plot([0.325 0.4 0.55 0.85],nocatch(2,:),'-','Color',map(6,:));
%     errbar([0.325 0.4 0.55 0.85],nocatch(2,:),nocatchs(2,:),'Color',map(6,:));


h1 = plot([0.15 0.3 0.45 0.6],nocatch(1,:),'o','MarkerSize',15);
set(h1(1),'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',map(1,:),'LineWidth',1.5);

h3 = plot(0.305,main(1),'o','MarkerSize',15);
set(h3(1),'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',map(2,:),'LineWidth',1.5);
%     errbar(0.305,main(1),mains(1),'Color',map(2,:));
h7 = plot(0.305,catch_(1),'o','MarkerSize',15);
%     errbar(0.305,catch_(1),catch_s(1),'Color',map(3,:));
set(h7(1),'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',map(3,:),'LineWidth',1.5);
h5 = plot([0.325 0.4 0.55 0.85],nocatch(2,:),'o','MarkerSize',15);
set(h5(1),'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',map(6,:),'LineWidth',1.5);
h9 = plot(0.4,main(2),'o','MarkerSize',15);
%     errbar(0.4,main(2),mains(2),'Color',map(5,:));
set(h9(1),'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',map(5,:),'LineWidth',1.5);

h4 = plot(0.4,catch_(2),'o','MarkerSize',15);
%     errbar(0.4,catch_(2),catch_s(2),'Color',map(4,:));
set(h4(1),'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',map(4,:),'LineWidth',1.5);


legend([h1,h3,h7,h5,h9,h4],{'Coherence: Control','Coherence: Attended','Coherence: Unattended','Contrast: Control','Contrast: Attended','Contrast: Unattended'});

title('Psychometric Functions for Cohcon');
xlabel('Contrast/Coherence (%)');
ylabel('Threshold (%)');
drawPublishAxis
%     set(h(1),'MarkerEdgeColor','r','MarkerFaceColor','none')
