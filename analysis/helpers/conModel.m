
function out = conModel(con,params,att,fitflag)

if isfield(params,'conmodel')
    if params.conmodel==1
        if isfield(params,'attgain') && params.attgain
            if att==1
                params.conslope = params.conslope * params.cohatt_congain;
            elseif att==2
                params.conslope = params.conslope * params.conatt_congain;
            end
        end
        out = params.conslope .* con;
    elseif params.conmodel==2
        if isfield(params,'attgain') && params.attgain
            if att==1
                params.conRmax = params.conRmax * params.cohatt_congain;
            elseif att==2
                params.conRmax = params.conRmax * params.conatt_congain;
            end
        end
        params.conn = round(params.conn);
        out = params.conRmax .* ((con.^params.conn) ./ (con.^params.conn + params.conc50.^params.conn));
    end
else
    warning('failure: no model selected');
    keyboard
end
if fitflag==0
    out = out + params.offset;
    if isfield(params,'attoff') && params.attoff
        if att==0
            return
        elseif att==1
            out = out + params.cohattoff;
        elseif att==2
            out = out + params.conattoff;
        else
            warning('failure');
            keyboard
        end
    end
end