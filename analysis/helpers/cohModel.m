
function out = cohModel(coh,params,att,fitflag)
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
    end
else
    warning('failure: no model selected');
    keyboard
end
if fitflag==0
    % setting fitflag = 0 means you're just trying to compute the response
    % model, but you're not fitting it. (This is so you can see the curves
    % and compare them)
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