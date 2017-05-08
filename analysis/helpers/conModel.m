function out = conModel(con,params)
global fixedParams

if isfield(params,'conmodel')
    if params.conmodel==1
        out = params.conslope .* con;
    elseif params.conmodel==2
        params.conn = round(params.conn);
        out = params.conRmax .* ((con.^params.conn) ./ (con.^params.conn + params.conc50.^params.conn));
    elseif params.conmodel==3 % exp model
        out = params.conalpha -(params.conalpha * exp(-params.conkappa*con));
    elseif params.conmodel==4
        % interpolation using existin gdata
        if isfield(params,'congain')
            lcon = fixedParams.con*params.congain;
        else
            lcon = fixedParams.con;
        end
        out = zeros(size(con));
        for ci = 1:length(con)
            out(ci) = lcon(find(fixedParams.x>=con(ci),1));
        end
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