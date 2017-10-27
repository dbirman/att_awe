function view_rca(out)
%% Startup
fname = out.filename;

% pull folder
cfolder = fname(1:(strfind(fname,'_data.mat')-1));

mrQuit
cd(fullfile('/Volumes/MoonData/data/cohcon_localizer/',cfolder));
folders = dir(pwd);

mrLoadRet

disp('When ready to continue type [dbcont]');
keyboard
%% Shutdown

mrQuit;