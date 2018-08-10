function out = conModel(con,params)
global fixedParams

if isfield(params,'deriv')
    if params.conmodel==2
        out = (params.conRmax * 1.9 * params.conc50^1.6 * con .^ 0.9 + 0.3 * con .^ 2.5) ./ ((params.conc50^params.conq+con.^params.conq).^2);
    end
    
    return
end

if isfield(params,'conmodel')
    if params.conmodel==1
        out = params.conslope .* con;
    elseif params.conmodel==2
        if isfield(params,'conn')
            out = params.conRmax .* ((con.^params.conn) ./ (con.^params.conn + params.conc50.^params.conn));
        else
            out = params.conRmax .* ((con.^(params.conp+params.conq)) ./ (con.^params.conq + params.conc50.^params.conq));
        end
    elseif params.conmodel==3 % exp model
        out = params.conalpha -(params.conalpha * exp(-params.conkappa*con));
    elseif params.conmodel==4
        % interpolation using existin gdata
        if isfield(params,'congain')
            lcon = fixedParams.con*params.congain;
        else
            lcon = fixedParams.con;
        end
        out = interp1(fixedParams.x,lcon,con);
        return
    elseif params.conmodel==5
        out = params.conalpha .* tanh(con*params.conkappa);
    elseif params.conmodel==6
        out = 1./exp(con).*(params.conalpha -(params.conalpha * exp(-params.conkappa*con))) + params.conalpha*con;
    end
else
    warning('failure: no model selected');
    keyboard
end