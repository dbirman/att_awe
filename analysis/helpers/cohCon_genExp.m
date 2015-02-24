function success = cohCon_genExp( inFile, outFile )
%GENEXP Converts stimulus matlab files to exp matlab files

try
    load(inFile);
    %% Peripheral re-org
    e = getTaskParameters(myscreen,task);
    expeye = getTaskEyeTraces(inFile,'taskNum=1','phaseNum=1','removeBlink=1');
    if isfield(expeye,'eye')
        e{2} = expeye.eye;
        clear expeye
    end
    
    %% Main Re-organization
    % for this experiment nothing...
    
    %% run Variables
    e{1}.runVars.pedestals = stimulus.pedestals;
    e{1}.runVars.runNum = stimulus.counter;
    e{1}.stimulus = stimulus;
    %% Save
    save(outFile,'e');
    success = 1;
catch
    success = 0;
end