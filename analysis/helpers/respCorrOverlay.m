function respCorrOverlay(out)
% RESPCORROVERLAY Generate the overlay for the analysis in out using the
% timeseries included.
%
%

%% Connect and load

fname = out.filename;

% pull folder
cfolder = fname(1:(strfind(fname,'_data.mat')-1));

mrQuit
cd(fullfile('/Volumes/MoonData/data/cohcon_localizer/',cfolder));
folders = dir(pwd);

%% Setup a view + Load Concatenation
view = newView();
view = viewSet(view,'curGroup','Concatenation');
view = viewSet(view,'curScan',1); % make sure scan # is correct

%% Load entire dataset
view = loadAnalysis(view,sprintf('erAnal/%s','all')); % check analysis name!

analysis = view.analyses{1};
tSeries = loadTSeries(view,view.curScan);

%% Correlation

tSeries_r = reshape(tSeries,size(tSeries,1)*size(tSeries,2)*size(tSeries,3),size(tSeries,4));
response = out.tSeries;

% next line crashes due to memory constraints
% r = corr(tSeries_r',response);

% use to avoid memory issues
disppercent(-1/size(tSeries_r,1));
for vi=1:size(tSeries_r,1)
    r(vi) = corr(out.tSeries,tSeries_r(vi,:)');
    disppercent(vi/size(tSeries_r,1));
end
disppercent(inf);

r_ = reshape(r,size(tSeries,1),size(tSeries,2),size(tSeries,3));

%% Save new overlays

mrDispOverlay({r_, r_.^2},view.curScan,view.curGroup,view,'overlayNames',{'r','r2'},'saveName=corrAnalysis_2');

%% Temp code: open mrLoadRet to view new overlays

mrQuit(0);