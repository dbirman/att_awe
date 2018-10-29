%%
%sids = {'s043','s300','s305','s329','s025'};
nSIDs = [305 329 43 25 300 346 343 344 338 340 348];
bSIDs = [345 330 337 335 349 354 353 356 334 352]; % behavior only participants
aSIDs = [nSIDs bSIDs];
sids = {};
for bi = 1:length(aSIDs)
    sids{end+1} = sprintf('s%03.0f',aSIDs(bi));
end
rois = {'V1','V2','V3','V4','V3a','V3b','V7','MT'};

%%
%%%%%%%%%%%%%%%%%%%%% NEURAL %%%%%%%%%%%%%%%%%%%%%%%%%
nSIDs = [305 329 43 25 300 346 343 338 344 348];
% NO ATT DATA FROM 340


%    1        2    3     4     5      6   7    8      9
% stimvol basecon lcon rcon basecoh lcoh rcoh timing task


%% Delete old files
% % warning('Are you sure?');
% % keyboard
% % for ni = 1:length(nSIDs)
% %     files = dir(fullfile(datafolder,sprintf('s%04.0f_a*.mat',nSIDs(ni))));
% %     for fi = 1:length(files)
% %         delete(fullfile(datafolder,files(fi).name));
% %     end
% % end
%%

for ni = 1:length(nSIDs)
    indiv_analysis(sprintf('s%04.0f',nSIDs(ni)),'decon');
end

%% Compare scanner behavior to outside scanner behavior


%% Cross analysis
% Cross analysis deconvolves the 16 combinations of contrast and coherence
% rather than separately deconvolving the contrast and coherence response.
% We rely on the rounding of the increments to do this, otherwise there's
% not really enough trials (<~25)
avg_deconResponseAtt_cross;
fit_deconResponseAtt_cross;
indiv_fitbehav_att_cross;

%% Sensitivity space to readout space analysis
sense_readout;