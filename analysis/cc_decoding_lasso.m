%%
fname = '~/proj/att_awe/analysis/data/s0300_voxels.csv';
[header,data] = csvreadh(fname);

%% Testing
aInst = aInsts{1}.l;

% get roi
roi_pos = find(cellfun(@(x) strcmp(x,inst.name),neural.ROIs));
croi = neural.ROIs{roi_pos};
[s, r] = parseROI(croi,neural.shortROIs);

% get stimvol
sN = neural.SCM_f.(folderz).lStim.stimNames;
sV = neural.SCM_f.(folderz).lStim.stimVol;
tsN = neural.SCM.(folderz).lStim.taskNames;
tsV = neural.SCM.(folderz).lStim.taskSV;
iV = inst.classify.instanceVol;
amps = inst.classify.instances;
[con,coh,task, sv] = sv2long(sV,sN,tsV,tsN);

% okay: organize like this:
%   Y          X1            X2
%  con     amp_vox_1     amp_vox_2
% etc...

% first, run through con, get the stimvol, find the matching amplitudes,
% and then expand

%%
tCount = 100000;
header = {'contrast','v1','v2','v3'};
data = zeros(tCount,106);
count = 1;

disp('Starting left...');
disppercent(-1/length(aInst));
for li = 1:length(aInst)
    inst = aInst{li};
    
    if ~isempty(inst.classify.instances)
        sN = neural.SCM_f.(folderz).lStim.stimNames;
        sV = neural.SCM_f.(folderz).lStim.stimVol;
        tsN = neural.SCM.(folderz).lStim.taskNames;
        tsV = neural.SCM.(folderz).lStim.taskSV;
        iV = inst.classify.instanceVol;
        amps = inst.classify.instances;
        [con,coh,task, sv] = sv2long(sV,sN,tsV,tsN);
        
        roi_pos = find(cellfun(@(x) strcmp(x,inst.name),neural.ROIs));
        croi = neural.ROIs{roi_pos};
        [s, r] = parseROI(croi,neural.shortROIs);
        
        
        for ci=1:length(amps)
            if ~isempty(amps)
                camps = amps{ci};
                for i = 1:size(camps,1)
                    camp = camps(i,:);
                    vol = iV{ci}(i);
                    
                    sv_pos = find(sv==vol);
                    if ~isempty(sv_pos)
                        ccon = con(sv_pos);
                        ccoh = coh(sv_pos);
                        ctask = task(sv_pos);
                        
                        dat = [ccon ccoh ctask fi s i camp];
                        
                        data(count,:) = dat;
                        count = count+1;
                    end
                end
                
            end
        end
    end
    
    data = data(1:count,:);
    
    disppercent(li/length(lInst));
end
disppercent(inf);

%% Lasso

Y = data(:,1); % contrast
X = data(:,7:80); % voxel amplitudes

Y = Y(~any(X>100,2),:);
X = X(~any(X>100,2),:); % remove rows with motion spikes first
X = X(:,~any(X>100,1)); % now remove voxels that show spikes

[B,stats] = lasso(X,Y,'CV',5);
lassoPlot(B,stats);

%% Get Betas

% B is voxels x lambda
Beta = B(:,stats.IndexMinMSE);
lassoPlot(B,stats,'PlotType','CV');
figure, hold on
plot(Y,Y,'*r') % X=Y plot)
plot(Y,X*Beta,'*b');
adjr2 = (std(Y)^2-stats.MSE(stats.IndexMinMSE))/(std(Y)^2);
disp(sprintf('Adjusted R^2: %0.2f',adjr2));

%%
Y2 = [X*Beta];
X2 = [ones(size(Y)),Y];
b2 = X2\Y2;
x = 0:.01:1;
plot(x,b2(1)+b2(2)*x,'*g');