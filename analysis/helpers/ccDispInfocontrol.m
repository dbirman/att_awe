
%%%%%%%%%%%%%%%%%%%%%%%
%    dispInfo    %
%%%%%%%%%%%%%%%%%%%%%%%
function ccDispInfocontrol(stimulus,subj)

%% No-Catch Performance
nocatch = zeros(2,4);
nocatchs = zeros(2,4);
params = zeros(2,4,4);

x = 0:.0001:1;

for task = 1:size(nocatch,1)
    for ped = 1:size(nocatch,2)
        try
        out = doStaircase('threshold',stimulus.staircases.nocatch{task,ped},'type','weibull','dispFig=0');
%         weib = weibull(x,out.fit.fitparams);
%         nocatch(task,ped) = x(find(weib>=normcdf(1/sqrt(2),0,1),1));
%             out = doStaircase('threshold',stimulus.staircases.nocatch{task,ped});
        nocatch(task,ped) = out.threshold;
        params(task,ped,:) = out.fit.fitparams;
%         nocatchs(task,ped) = out.thresholdSTE;
        catch % probably a missing staircase, whatever
        end
    end
end
if any(nocatch(:)>10)
    warning('removing stairs >10');
    nocatch(nocatch>10) = 0;
end

%% Save data
fname = fullfile(datafolder,sprintf('%s_data.mat',subj));
if exist(fname)==2, load(fname); end
data.control = nocatch;
data.params = params;
save(fname,'data');