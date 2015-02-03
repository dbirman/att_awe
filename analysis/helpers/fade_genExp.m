function success = fade_genExp( inFile, outFile )
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
    e{1}.randVars.image1 = zeros(size(e{1}.randVars.deltaPed));
    e{1}.randVars.image2 = zeros(size(e{1}.randVars.deltaPed));
    e{1}.randVars.image3 = zeros(size(e{1}.randVars.deltaPed));
    e{1}.randVars.image4 = zeros(size(e{1}.randVars.deltaPed));
    e{1}.randVars.gender1 = zeros(size(e{1}.randVars.deltaPed));
    e{1}.randVars.gender2 = zeros(size(e{1}.randVars.deltaPed));
    e{1}.randVars.gender3 = zeros(size(e{1}.randVars.deltaPed));
    e{1}.randVars.gender4 = zeros(size(e{1}.randVars.deltaPed));
    e{1}.randVars.contrast1 = zeros(size(e{1}.randVars.deltaPed));
    e{1}.randVars.contrast2 = zeros(size(e{1}.randVars.deltaPed));
    e{1}.randVars.contrast3 = zeros(size(e{1}.randVars.deltaPed));
    e{1}.randVars.contrast4 = zeros(size(e{1}.randVars.deltaPed));
    e{1}.randVars.tNoise = zeros(size(e{1}.randVars.deltaPed));
    e{1}.randVars.tCon = zeros(size(e{1}.randVars.deltaPed));
    e{1}.randVars.tGen = zeros(size(e{1}.randVars.deltaPed));
    e{1}.randVars.tImage = zeros(size(e{1}.randVars.deltaPed));
    for i = 1:length(e{1}.randVars.deltaPed)
        cuePos = e{1}.randVars.target(i);
        changePos = e{1}.randVars.changeTarget(i);
%         e{1}.randVars.tNoise(i) = e{1}.randVars.noiseList{i}(cuePos,1);
        e{1}.randVars.tCon(i) = e{1}.randVars.contrastList{i}(cuePos);
%         e{1}.randVars.tGen(i) = e{1}.randVars.genderList{i}(cuePos);
        e{1}.randVars.tImage(i) = e{1}.randVars.imageNums{i}(cuePos);
%         e{1}.randVars.cNoise(i) = e{1}.randVars.noiseList{i}(changePos,1);
        e{1}.randVars.cCon(i) = e{1}.randVars.contrastList{i}(changePos);
%         e{1}.randVars.cGen(i) = e{1}.randVars.genderList{i}(changePos);
        e{1}.randVars.cImage(i) = e{1}.randVars.imageNums{i}(changePos);
        e{1}.randVars.image1(i) = e{1}.randVars.imageNums{i}(1);
        e{1}.randVars.image2(i) = e{1}.randVars.imageNums{i}(2);
        e{1}.randVars.image3(i) = e{1}.randVars.imageNums{i}(3);
        e{1}.randVars.image4(i) = e{1}.randVars.imageNums{i}(4);
%         e{1}.randVars.gender1(i) = e{1}.randVars.genderList{i}(1);
%         e{1}.randVars.gender2(i) = e{1}.randVars.genderList{i}(2);
%         e{1}.randVars.gender3(i) = e{1}.randVars.genderList{i}(3);
%         e{1}.randVars.gender4(i) = e{1}.randVars.genderList{i}(4);
        e{1}.randVars.contrast1 = e{1}.randVars.contrastList{i}(1,1);
        e{1}.randVars.contrast2 = e{1}.randVars.contrastList{i}(1,1);
        e{1}.randVars.contrast3 = e{1}.randVars.contrastList{i}(1,1);
        e{1}.randVars.contrast4 = e{1}.randVars.contrastList{i}(1,1);
% %         e{1}.randVars.noise1_int1(i) = e{1}.randVars.noiseList{i}(1,1);
% %         e{1}.randVars.noise1_int2(i) = e{1}.randVars.noiseList{i}(1,2);
% %         e{1}.randVars.noise2_int1(i) = e{1}.randVars.noiseList{i}(2,1);
% %         e{1}.randVars.noise2_int2(i) = e{1}.randVars.noiseList{i}(2,2);
    end

    %% run Variables
    e{1}.runVars.contrastLists = stimulus.blocks.contrastListOptions(stimulus.blocks.counter-6:stimulus.blocks.counter,:);
    e{1}.runVars.pedestals = stimulus.pedestals;
    e{1}.runVars.runNum = stimulus.counter;
    e{1}.stimulus = stimulus;
    %% Save
    save(outFile,'e');
    success = 1;
catch
    success = 0;
end