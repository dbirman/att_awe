function fit = fit_doublegamma(hrf,time)

global params fixedParams

hrfparams.amp1 = 1;
hrfparams.tau1 = [0.4 -inf inf];
hrfparams.timelag1 = [1.5 0 3];
hrfparams.amp2 = [-0.2 -inf 0];
hrfparams.tau2 = [0.2];
hrfparams.timelag2 = [4];
hrfparams.exponent = 7;
hrfparams.offset = [0 0 inf];

params.hrfparams = hrfparams;
params.roiparams 

[initparams,minparams,maxparams] = initParams();

optimParams = optimset('Algorithm','trust-region-reflective','MaxIter',inf,'Display','off');

[bestparams, ~, ~, ~, ~, ~, ~] = lsqnonlin(@hrfResidual,initparams,minparams,maxparams,optimParams,time,hrf);

fit.params = getParams(bestparams,fixedParams);
fit.time = time;
fit.hrf = hrf;
[res,out] = hrfResidual(bestparams,time,hrf);
fit.out = out;
fit.res = res;

function [res,out] = hrfResidual(params,time,hrf)

global fixedParams

params = getParams(params,fixedParams);

out = cc_gamma(time,params);

res = out-hrf;


function [initparams, minparams, maxparams] = initParams()
%%
global fixedParams params

%% Deal with HRF params
fixedParams.strs = fields(params.hrfparams)';

initparams = [];
minparams = [];
maxparams = [];
indexes = {};
count = 1;

fixed = [];
optim = zeros(size(indexes));

count = 1;
for i = 1:length(fixedParams.strs)
    cvals = params.hrfparams.(fixedParams.strs{i});
    if length(cvals)==1
        fixedParams.(fixedParams.strs{i}) = cvals;
        fixed(end+1) = 1;
        indexes{end+1} = [];
    elseif length(cvals)==3
        initparams(end+1) = cvals(1);
        minparams(end+1) = cvals(2);
        maxparams(end+1) = cvals(3);
        fixed(end+1) = 0;
        indexes{end+1} = count; count = count+1;
    elseif length(cvals)==2 || length(cvals)>3
        warning('Failure');
        keyboard
        % optimizer
        fixedParams.(fixedParams.strs{i}) = cvals;
        optim(i) = 1;
    else
        error('You initialized a parameter with the wrong initial values... unable to interpret');
    end
end
% 
% %% Deal with ROI params
% if iscell(params.roiparams)
%     rStrs = {};
%     disp('ROI parameters were passed in as a cell: interpreting as fixed values');
%     for ri = 1:length(fixedParams.ROIs)
%         % get the fields
%         crfields = fields(params.roiparams{ri});
%         % copy each field with this ROIs prefix
%         for ci = 1:length(crfields)
%             rStrs{end+1} = sprintf('%s%s',fixedParams.ROIs{ri},crfields{ci});
%             fixedParams.(rStrs{end}) = params.roiparams{ri}.(crfields{ci});
%             fixed(end+1) = 1;
%             indexes{end+1} = [];
%         end
%     end
% else
%     rfields = fields(params.roiparams);
% 
%     rStrs = {};
%     for ni = 1:length(rfields)
%         cvals = params.roiparams.(rfields{ni});
%         % check if this name includes an ROI already, if so, just copy it in
%         % directly
%         cfield = rfields{ni};
%         if any(cellfun(@(x) ~isempty(strfind(cfield,x)),fixedParams.ROIs))
%             % just copy directly
%             rStrs{end+1} = cfield;
%             fixedParams.(rStrs{end}) = cvals;
%             fixed(end+1) = 1;
%             indexes{end+1} = [];
%         else
%             % replicate for every ROI
%             for ri = 1:length(fixedParams.ROIs)
%                 rStrs{end+1} = sprintf('%s%s',fixedParams.ROIs{ri},rfields{ni});
%                 if length(cvals)==1
%                     fixedParams.(rStrs{end}) = cvals;
%                     fixed(end+1) = 1;
%                     indexes{end+1} = [];
%                 elseif length(cvals)==3
%                     initparams(end+1) = cvals(1);
%                     minparams(end+1) = cvals(2);
%                     maxparams(end+1) = cvals(3);
%                     indexes{end+1} = count;
%                     count = count+1;
%                     fixed(end+1) = 0;
%                 elseif length(cvals)==2 || length(cvals)>3
%                     error('You are not allowed to use the optimizer for ROI specific parameters...');
%                 else
%                     error('You initialized a parameter with the wrong initial values... unable to interpret');
%                 end
%             end
%         end
%     end 
% end

%% Save optim/fixed/indexes
% fixedParams.strs = [fixedParams.strs rStrs];
fixedParams.optim = optim;
fixedParams.fixed = fixed;
fixedParams.idx = indexes;
% 
% function p = getROIParams(params,ROI)
% % just grab anything that starts with ROI
% p = struct;
% flds = fields(params);
% for fi = 1:length(flds)
%     if strfind(flds{fi},ROI)
%         p.(strrep(flds{fi},ROI,'')) = params.(flds{fi});
%     end
% end

function p = getParams(params,fixedParams)

for i = 1:length(fixedParams.strs)
    if fixedParams.fixed(i)
        p.(fixedParams.strs{i}) = fixedParams.(fixedParams.strs{i});
    else
        p.(fixedParams.strs{i}) = params(fixedParams.idx{i});
    end
end