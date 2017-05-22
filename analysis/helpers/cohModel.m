function out = cohModel(coh,params)
global fixedParams

if isfield(params,'cohmodel')
    if params.cohmodel==1
        out = params.cohslope .* coh;
    elseif params.cohmodel==2
%         out = params.cohRmax .* ((coh.^params.cohn) ./ (coh.^params.cohn + params.cohc50.^params.cohn));
        out = params.cohRmax .* ((coh.^(params.cohp+params.cohq)) ./ (coh.^params.cohq + params.cohc50.^params.cohq));
    elseif params.cohmodel ==3
        out = params.cohalpha-(params.cohalpha * exp(-params.cohkappa*coh));
    elseif params.cohmodel==4
        if isfield(params,'cohgain')
            lcoh = fixedParams.coh*params.cohgain;
        else
            lcoh = fixedParams.coh;
        end
        out = zeros(size(coh));
        for ci = 1:length(coh)
            out(ci) = lcoh(find(fixedParams.x>=coh(ci),1));
        end
        return
    elseif params.cohmodel==5
        out = params.cohalpha .* tanh(coh*params.cohkappa);
    elseif params.cohmodel==6
        out = 1./exp(coh).*(params.cohalpha -(params.cohalpha * exp(-params.cohkappa*coh)))+params.cohalpha*coh;
    end
else
    warning('failure: no model selected');
    keyboard
end