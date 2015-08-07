function SCM = cc_removeNoDataStimVols(SCM,name,thresh,restore)

if restore
    data = SCM.(name);

    scmFields = fields(data);

    for fi = 1:length(scmFields)
        cfield = scmFields{fi};

        data.(cfield) = data.(cfield).backup;
    end

    SCM.(name) = data;
    return
end

data = SCM.(name);

scmFields = fields(data);

for fi = 1:length(scmFields)
    cfield = scmFields{fi};
    
    data.(cfield).backup = data.(cfield);
    
    dat = data.(cfield);
    
    i = 1;
    while i < length(dat.stimVol)
        if length(dat.stimVol{i}) < thresh
            disp(sprintf('(removeNoSV) Removing %i volumes for stimulus: %s',length(dat.stimVol{i}),dat.stimNames{i}));
            dat.stimNames(i) = [];
            dat.stimVol(i) = [];
        else
            i = i + 1;
        end
    end
    
    data.(cfield) = dat;
end

SCM.(name) = data;