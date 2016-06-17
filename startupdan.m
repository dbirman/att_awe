dbstop('if','error');

pathNames = {'~/proj/att_awe'};

for i = 1:length(pathNames)
  if isdir(pathNames{i})
    addpath(genpath_exclude(pathNames{i},{'.git','.svm'}));
  end
end

mrSetPref('maxBlocksize',32000000000);
mglSetParam('sidDatabaseFilename','/Volumes/gru/sid/sidDatabase')

% run user specific matlab startup. This will look for a matlab
% script like startupkenji to run, if it finds it on the path
% will run
% username = getusername;
% userStartup = sprintf('startup%s',username);
% if exist(userStartup) == 2
%   disp(sprintf('(startup) Running user specific startup: %s',userStartup));
%   eval(userStartup);
% else
%   disp(sprintf('(startup) No user specific startup scrip found. If you need to, create one called: %s',userStartup));
% end  

cd('~/proj/')