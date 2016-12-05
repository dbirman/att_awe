function out = conModel(con,params,att,fitflag)
global fixedParams

if ~exist('att','var')
    att=0;
end
if ~exist('fitflag','var')
    fitflag=0;
end
if isfield(params,'attgain')
    warning('Attention gain model failure');
end

if isfield(params,'conmodel')
    if params.conmodel==1
%         if isfield(params,'attgain') && params.attgain
%             if att==1
%                 params.conslope = params.conslope * params.cohatt_congain;
%             elseif att==2
%                 params.conslope = params.conslope * params.conatt_congain;
%             end
%         end
        out = params.conslope .* con;
    elseif params.conmodel==2
%         if isfield(params,'attgain') && params.attgain
%             if att==1
%                 params.conRmax = params.conRmax * params.cohatt_congain;
%             elseif att==2
%                 params.conRmax = params.conRmax * params.conatt_congain;
%             end
%         end
        params.conn = round(params.conn);
        out = params.conRmax .* ((con.^params.conn) ./ (con.^params.conn + params.conc50.^params.conn));
    elseif params.conmodel==3 % exp model
        out = -params.conalpha * exp(-params.conkappa*con);
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
    end
else
    warning('failure: no model selected');
    keyboard
end

if isfield(params,'conalpha')
    out = out+params.conalpha; % add the alpha parameter so that the function starts at zero
end

return
% if ~isfield(fixedParams,'fitting') || fixedParams.fitting==0
%     out = out + params.offset; % adds the offset, it'll do this equally for
%     % contrast and coherence which isn't strictly true, but for
%     % visualization purposes it's a reasonable approximation
%     
%     if isfield(params,'attoff') && params.attoff
%         warning('code not implemented!!');
%         keyboard
%         if att==0
%             return
%         elseif att==1
%             out = out + params.cohattoff;
%         elseif att==2
%             out = out + params.conattoff;
%         else
%             warning('failure');
%             keyboard
%         end
%     end
% end