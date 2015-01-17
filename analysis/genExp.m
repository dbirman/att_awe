function success = genExp( inFile, outFile )
%GENEXP Converts stimulus matlab files to exp matlab files

try
    load(inFile);
    %% Peripheral re-org
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
    %% Main Re-organization
    exp{1}.correct = exp{1}.response==exp{1}.randVars.interval;
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
    exp{1}.randVars.tNoise1 = zeros(size(exp{1}.randVars.deltaPed));
    exp{1}.randVars.tNoise2 = zeros(size(exp{1}.randVars.deltaPed));
    exp{1}.randVars.tCon1 = zeros(size(exp{1}.randVars.deltaPed));
    exp{1}.randVars.tCon2 = zeros(size(exp{1}.randVars.deltaPed));
    exp{1}.randVars.tGen = zeros(size(exp{1}.randVars.deltaPed));
    exp{1}.randVars.tImage = zeros(size(exp{1}.randVars.deltaPed));
    for i = 1:length(exp{1}.randVars.deltaPed)
        pos = exp{1}.randVars.targetLoc(i);
        exp{1}.randVars.tNoise1 = exp{1}.randVars.noiseList{i}(pos,1);
        exp{1}.randVars.tNoise2 = exp{1}.randVars.noiseList{i}(pos,2);
        exp{1}.randVars.tCon1 = exp{1}.randVars.contrastList{i}(pos,1);
        exp{1}.randVars.tCon2 = exp{1}.randVars.contrastList{i}(pos,2);
        exp{1}.randVars.tGen = exp{1}.randVars.genderList{i}(pos);
        exp{1}.randVars.tImage = exp{1}.randVars.imageNums{i}(pos);
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
    %% Staircase Fix:
    % Staircase is, (NOISE1/CON2),CUE1/CUE4,PED1:3
    exp{1}.stairNoise = squeeze(stimulus.staircase(1,:,:));
    exp{1}.stairNoiseDual = squeeze(stimulus.dualstaircase(1,:,:));
    exp{1}.stairCon = squeeze(stimulus.staircase(2,:,:));
    exp{1}.stairConDual = squeeze(stimulus.dualstaircase(2,:,:));
    exp{2}.stair = squeeze(stimulus.p.staircase);
    exp{2}.stairDual = squeeze(stimulus.p.dualstaircase);
    %% Peripheral Re-organization
    %           f m non
    per_stim = [9 8 10];
    exp{2}.randVars.tGen = zeros(size(exp{2}.randVars.position));
    exp{2}.randVars.tImage = zeros(size(exp{2}.randVars.position));
    exp{2}.correct = zeros(size(exp{2}.randVars.position));
    for i = 1:length(exp{2}.randVars.position)
        pos = exp{2}.randVars.position(i);
        exp{2}.randVars.tGen(i) = exp{2}.randVars.gender{i}(pos);
        exp{2}.randVars.tImage(i) = exp{2}.randVars.images{i}(pos);
        if isnan(exp{2}.response(i))
            exp{2}.correct(i) = nan;
        else
            if exp{2}.response(i) == 10 && exp{2}.randVars.respond(i) == 0
                exp{2}.correct(i) = 1;
            elseif exp{2}.randVars.respond(i) == 1
                exp{2}.correct(i) = exp{2}.response(i)==per_stim(exp{2}.randVars.tGen(i));
            else
                exp{2}.correct(i) = 0;
            end
        end
    end
    exp{2}.randVars = rmfield(exp{2}.randVars,'gender');
    exp{2}.randVars = rmfield(exp{2}.randVars,'images');
    %% run Variables
    exp{3}.runVars.blocks = stimulus.blockList;
    exp{3}.runVars.dual = stimulus.dualList(end);
    exp{3}.runVars.pedestals = stimulus.pedestals;
    exp{3}.runVars.runNum = stimulus.counter;
    exp{3}.stimulus = stimulus;
    %% Save
    save(outFile,'exp');
    success = 1;
catch
    success = 0;
end