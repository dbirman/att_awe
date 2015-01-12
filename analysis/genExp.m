function success = genExp( inFile, outFile )
%GENEXP Converts stimulus matlab files to exp matlab files

try
    load(inFile);
    exp = getTaskParameters(myscreen,task);
    expeye = getTaskEyeTraces(inFile,'taskNum=1','phaseNum=1','segNum=2','removeBlink=1');
    pereye = getTaskEyeTraces(inFile,'taskNum=2','phaseNum=1','segNum=2','removeBlink=1');
    if isfield(expeye,'eye')
        exp.eye = expeye.eye;
        clear expeye
    end
    if isfield(pereye,'eye')
        exp.p.eye = pereye.eye;
        clear pereye;
    end
    %% Peripheral Re-organization
    exp{2}.randVars.gender1 = zeros(size(exp{2}.randVars.position));
    exp{2}.randVars.gender2 = zeros(size(exp{2}.randVars.position));
    exp{2}.randVars.image1 = zeros(size(exp{2}.randVars.position));
    exp{2}.randVars.image2 = zeros(size(exp{2}.randVars.position));
    for i = 1:length(exp{2}.randVars.position)
        exp{2}.randVars.gender1(i) = exp{2}.randVars.gender{i}(1);
        exp{2}.randVars.gender2(i) = exp{2}.randVars.gender{i}(2);
        exp{2}.randVars.image1(i) = exp{2}.randVars.images{i}(2);
        exp{2}.randVars.image2(i) = exp{2}.randVars.images{i}(2);
    end
    exp{2}.randVars = rmfield(exp{2}.randVars,'gender');
    exp{2}.randVars = rmfield(exp{2}.randVars,'images');
    %% Main Re-organization
    exp{1}.randVars.image1 = zeros(size(exp{1}.randVars.deltaPed));
    exp{1}.randVars.image2 = zeros(size(exp{1}.randVars.deltaPed));
    exp{1}.randVars.image3 = zeros(size(exp{1}.randVars.deltaPed));
    exp{1}.randVars.image4 = zeros(size(exp{1}.randVars.deltaPed));
    exp{1}.randVars.gender1 = zeros(size(exp{1}.randVars.deltaPed));
    exp{1}.randVars.gender2 = zeros(size(exp{1}.randVars.deltaPed));
    exp{1}.randVars.gender3 = zeros(size(exp{1}.randVars.deltaPed));
    exp{1}.randVars.gender4 = zeros(size(exp{1}.randVars.deltaPed));
    exp{1}.randVars.contrast1_int1 = zeros(size(exp{1}.randVars.deltaPed));
    exp{1}.randVars.contrast1_int2 = zeros(size(exp{1}.randVars.deltaPed));
    exp{1}.randVars.contrast2_int1 = zeros(size(exp{1}.randVars.deltaPed));
    exp{1}.randVars.contrast2_int2 = zeros(size(exp{1}.randVars.deltaPed));
    exp{1}.randVars.contrast3_int1 = zeros(size(exp{1}.randVars.deltaPed));
    exp{1}.randVars.contrast3_int2 = zeros(size(exp{1}.randVars.deltaPed));
    exp{1}.randVars.contrast4_int1 = zeros(size(exp{1}.randVars.deltaPed));
    exp{1}.randVars.contrast4_int2 = zeros(size(exp{1}.randVars.deltaPed));
    exp{1}.randVars.noise1_int1 = zeros(size(exp{1}.randVars.deltaPed));
    exp{1}.randVars.noise1_int2 = zeros(size(exp{1}.randVars.deltaPed));
    exp{1}.randVars.noise2_int1 = zeros(size(exp{1}.randVars.deltaPed));
    exp{1}.randVars.noise2_int2 = zeros(size(exp{1}.randVars.deltaPed));
    exp{1}.randVars.noise3_int1 = zeros(size(exp{1}.randVars.deltaPed));
    exp{1}.randVars.noise3_int2 = zeros(size(exp{1}.randVars.deltaPed));
    exp{1}.randVars.noise4_int1 = zeros(size(exp{1}.randVars.deltaPed));
    exp{1}.randVars.noise4_int2 = zeros(size(exp{1}.randVars.deltaPed));
    for i = 1:length(exp{1}.randVars.deltaPed)
        exp{1}.randVars.image1(i) = exp{1}.randVars.imageNums{i}(1);
        exp{1}.randVars.image2(i) = exp{1}.randVars.imageNums{i}(2);
        exp{1}.randVars.image3(i) = exp{1}.randVars.imageNums{i}(3);
        exp{1}.randVars.image4(i) = exp{1}.randVars.imageNums{i}(4);
        exp{1}.randVars.gender1(i) = exp{1}.randVars.genderList{i}(1);
        exp{1}.randVars.gender2(i) = exp{1}.randVars.genderList{i}(2);
        exp{1}.randVars.gender3(i) = exp{1}.randVars.genderList{i}(3);
        exp{1}.randVars.gender4(i) = exp{1}.randVars.genderList{i}(4);
        exp{1}.randVars.contrast1_int1(i) = exp{1}.randVars.contrastList{i}(1,1);
        exp{1}.randVars.contrast1_int2(i) = exp{1}.randVars.contrastList{i}(1,2);
        exp{1}.randVars.contrast2_int1(i) = exp{1}.randVars.contrastList{i}(2,1);
        exp{1}.randVars.contrast2_int2(i) = exp{1}.randVars.contrastList{i}(2,2);
        exp{1}.randVars.contrast3_int1(i) = exp{1}.randVars.contrastList{i}(3,1);
        exp{1}.randVars.contrast3_int2(i) = exp{1}.randVars.contrastList{i}(3,2);
        exp{1}.randVars.contrast4_int1(i) = exp{1}.randVars.contrastList{i}(4,1);
        exp{1}.randVars.contrast4_int2(i) = exp{1}.randVars.contrastList{i}(4,2);
        exp{1}.randVars.noise1_int1(i) = exp{1}.randVars.noiseList{i}(1,1);
        exp{1}.randVars.noise1_int2(i) = exp{1}.randVars.noiseList{i}(1,2);
        exp{1}.randVars.noise2_int1(i) = exp{1}.randVars.noiseList{i}(2,1);
        exp{1}.randVars.noise2_int2(i) = exp{1}.randVars.noiseList{i}(2,2);
        exp{1}.randVars.noise3_int1(i) = exp{1}.randVars.noiseList{i}(3,1);
        exp{1}.randVars.noise3_int2(i) = exp{1}.randVars.noiseList{i}(3,2);
        exp{1}.randVars.noise4_int1(i) = exp{1}.randVars.noiseList{i}(4,1);
        exp{1}.randVars.noise4_int2(i) = exp{1}.randVars.noiseList{i}(4,2);
    end
    exp{1}.randVars = rmfield(exp{1}.randVars,'imageNums');
    exp{1}.randVars = rmfield(exp{1}.randVars,'genderList');
    exp{1}.randVars = rmfield(exp{1}.randVars,'contrastList');
    exp{1}.randVars = rmfield(exp{1}.randVars,'noiseList');
    %% Stimulus Fix:
    % Staircase is, (NOISE1/CON2),CUE1/CUE4,PED1:3
    exp{1}.stairNoise = squeeze(stimulus.staircase(1,:,:));
    exp{1}.stairNoiseDual = squeeze(stimulus.dualstaircase(1,:,:));
    exp{1}.stairCon = squeeze(stimulus.staircase(2,:,:));
    exp{1}.stairConDual = squeeze(stimulus.dualstaircase(2,:,:));
    exp{1}.p.stair = squeeze(stimulus.p.staircase);
    exp{1}.p.stairDual = squeeze(stimulus.p.dualstaircase);
    %% run Variables
    exp{3}.runVars.blocks = stimulus.blockList;
    exp{3}.runVars.dual = stimulus.dualList;
    exp{3}.runVars.pedestals = stimulus.pedestals;
    %% Save
    save(outFile,'exp');
    success = 1;
catch
    success = 0;
end