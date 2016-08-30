
function out = cohModel(coh,params,att,fitflag)
global fixedParams

if ~exist('att','var')
    att=0;
end
if ~exist('fitflag','var')
    fitflag=0;
end
if isfield(params,'cohmodel')
    if params.cohmodel==1
        if isfield(params,'attgain') && params.attgain
            if att==1
                params.cohslope = params.cohslope * params.cohatt_cohgain;
            elseif att==2
                params.cohslope = params.cohslope * params.conatt_cohgain;
            end
        end
        out = params.cohslope .* coh;
    elseif params.cohmodel==2
        if isfield(params,'attgain') && params.attgain==1
            if att==1
                params.cohRmax = params.cohRmax * params.cohatt_cohgain;
            elseif att==2
                params.cohRmax = params.cohRmax * params.conatt_cohgain;
            end
        end
        params.cohn = round(params.cohn);
        out = params.cohRmax .* ((coh.^params.cohn) ./ (coh.^params.cohn + params.cohc50.^params.cohn));
    elseif params.cohmodel ==3
        if isfield(params,'attgain')
            warning('attgain not implemented for exponential model');
        end
        out = -params.cohalpha * exp(-params.cohkappa*coh);
    elseif params.cohmodel==4
        out = zeros(size(coh));
        for ci = 1:length(coh)
            out(ci) = fixedParams.coh(find(fixedParams.x>=coh(ci),1));
        end
        return
    end
else
    warning('failure: no model selected');
    keyboard
end
if fitflag==0
    % setting fitflag = 0 means you're just trying to compute the response
    % model, but you're not fitting it. (This is so you can see the curves
    % and compare them)
    if ~isfield(params,'offset')
        % this is probably the behavioral model, just return
        if isfield(params,'cohalpha')
            out = out + params.cohalpha;
        end
        return
    end
    out = out + params.offset;
    if isfield(params,'attoff') && params.attoff
        if att==0
            return
        elseif att==2
            out = out + params.conattoff;
        elseif att==1
            out = out + params.cohattoff;
        else
            warning('failure');
            keyboard
        end
    end
end