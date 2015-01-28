function success = genExp( inFile, outFile )
%GENEXP Converts stimulus matlab files to exp matlab files

try
    load(inFile);
    %% Peripheral re-org
    e = getTaskParameters(myscreen,task);
    expeye = getTaskEyeTraces(inFile,'taskNum=1','phaseNum=1','removeBlink=1');
    pereye = getTaskEyeTraces(inFile,'taskNum=2','phaseNum=1','segNum=2','removeBlink=1');
    if isfield(expeye,'eye')
        e{4} = expeye.eye;
        clear expeye
    end
    if isfield(pereye,'eye')
        e{5} = pereye.eye;
        clear pereye;
    end
    
    %% Add timing variables to peripheral
    for t = 1:length(e{2}.trials)
        if e{2}.randVars.respond(t)
            time = e{2}.trials(t).segtime;
            try
                e{2}.stimPresLength(t) = time(4)-time(3);
            catch
                e{2}.stimPresLength(t) = NaN;
            end
        else
            e{2}.stimPresLength(t) = NaN;
        end
    end
    
    %% Main Re-organization
    e{1}.correct = e{1}.response==e{1}.randVars.interval;
    e{1}.randVars.image1 = zeros(size(e{1}.randVars.deltaPed));
    e{1}.randVars.image2 = zeros(size(e{1}.randVars.deltaPed));
    e{1}.randVars.image3 = zeros(size(e{1}.randVars.deltaPed));
    e{1}.randVars.image4 = zeros(size(e{1}.randVars.deltaPed));
    e{1}.randVars.gender1 = zeros(size(e{1}.randVars.deltaPed));
    e{1}.randVars.gender2 = zeros(size(e{1}.randVars.deltaPed));
    e{1}.randVars.gender3 = zeros(size(e{1}.randVars.deltaPed));
    e{1}.randVars.gender4 = zeros(size(e{1}.randVars.deltaPed));
    e{1}.randVars.contrast1_int1 = zeros(size(e{1}.randVars.deltaPed));
    e{1}.randVars.contrast1_int2 = zeros(size(e{1}.randVars.deltaPed));
    e{1}.randVars.contrast2_int1 = zeros(size(e{1}.randVars.deltaPed));
    e{1}.randVars.contrast2_int2 = zeros(size(e{1}.randVars.deltaPed));
    e{1}.randVars.contrast3_int1 = zeros(size(e{1}.randVars.deltaPed));
    e{1}.randVars.contrast3_int2 = zeros(size(e{1}.randVars.deltaPed));
    e{1}.randVars.contrast4_int1 = zeros(size(e{1}.randVars.deltaPed));
    e{1}.randVars.contrast4_int2 = zeros(size(e{1}.randVars.deltaPed));
    e{1}.randVars.noise1_int1 = zeros(size(e{1}.randVars.deltaPed));
    e{1}.randVars.noise1_int2 = zeros(size(e{1}.randVars.deltaPed));
    e{1}.randVars.noise2_int1 = zeros(size(e{1}.randVars.deltaPed));
    e{1}.randVars.noise2_int2 = zeros(size(e{1}.randVars.deltaPed));
    e{1}.randVars.noise3_int1 = zeros(size(e{1}.randVars.deltaPed));
    e{1}.randVars.noise3_int2 = zeros(size(e{1}.randVars.deltaPed));
    e{1}.randVars.noise4_int1 = zeros(size(e{1}.randVars.deltaPed));
    e{1}.randVars.noise4_int2 = zeros(size(e{1}.randVars.deltaPed));
    e{1}.randVars.tNoise1 = zeros(size(e{1}.randVars.deltaPed));
    e{1}.randVars.tNoise2 = zeros(size(e{1}.randVars.deltaPed));
    e{1}.randVars.tCon1 = zeros(size(e{1}.randVars.deltaPed));
    e{1}.randVars.tCon2 = zeros(size(e{1}.randVars.deltaPed));
    e{1}.randVars.tGen = zeros(size(e{1}.randVars.deltaPed));
    e{1}.randVars.tImage = zeros(size(e{1}.randVars.deltaPed));
    for i = 1:length(e{1}.randVars.deltaPed)
        pos = e{1}.randVars.targetLoc(i);
        e{1}.randVars.tNoise1(i) = e{1}.randVars.noiseList{i}(pos,1);
        e{1}.randVars.tNoise2(i) = e{1}.randVars.noiseList{i}(pos,2);
        e{1}.randVars.tCon1(i) = e{1}.randVars.contrastList{i}(pos,1);
        e{1}.randVars.tCon2(i) = e{1}.randVars.contrastList{i}(pos,2);
        e{1}.randVars.tGen(i) = e{1}.randVars.genderList{i}(pos);
        e{1}.randVars.tImage(i) = e{1}.randVars.imageNums{i}(pos);
        e{1}.randVars.image1(i) = e{1}.randVars.imageNums{i}(1);
        e{1}.randVars.image2(i) = e{1}.randVars.imageNums{i}(2);
        e{1}.randVars.image3(i) = e{1}.randVars.imageNums{i}(3);
        e{1}.randVars.image4(i) = e{1}.randVars.imageNums{i}(4);
        e{1}.randVars.gender1(i) = e{1}.randVars.genderList{i}(1);
        e{1}.randVars.gender2(i) = e{1}.randVars.genderList{i}(2);
        e{1}.randVars.gender3(i) = e{1}.randVars.genderList{i}(3);
        e{1}.randVars.gender4(i) = e{1}.randVars.genderList{i}(4);
        e{1}.randVars.contrast1_int1(i) = e{1}.randVars.contrastList{i}(1,1);
        e{1}.randVars.contrast1_int2(i) = e{1}.randVars.contrastList{i}(1,2);
        e{1}.randVars.contrast2_int1(i) = e{1}.randVars.contrastList{i}(2,1);
        e{1}.randVars.contrast2_int2(i) = e{1}.randVars.contrastList{i}(2,2);
        e{1}.randVars.contrast3_int1(i) = e{1}.randVars.contrastList{i}(3,1);
        e{1}.randVars.contrast3_int2(i) = e{1}.randVars.contrastList{i}(3,2);
        e{1}.randVars.contrast4_int1(i) = e{1}.randVars.contrastList{i}(4,1);
        e{1}.randVars.contrast4_int2(i) = e{1}.randVars.contrastList{i}(4,2);
        e{1}.randVars.noise1_int1(i) = e{1}.randVars.noiseList{i}(1,1);
        e{1}.randVars.noise1_int2(i) = e{1}.randVars.noiseList{i}(1,2);
        e{1}.randVars.noise2_int1(i) = e{1}.randVars.noiseList{i}(2,1);
        e{1}.randVars.noise2_int2(i) = e{1}.randVars.noiseList{i}(2,2);
        e{1}.randVars.noise3_int1(i) = e{1}.randVars.noiseList{i}(3,1);
        e{1}.randVars.noise3_int2(i) = e{1}.randVars.noiseList{i}(3,2);
        e{1}.randVars.noise4_int1(i) = e{1}.randVars.noiseList{i}(4,1);
        e{1}.randVars.noise4_int2(i) = e{1}.randVars.noiseList{i}(4,2);
    end
    e{1}.randVars = rmfield(e{1}.randVars,'imageNums');
    e{1}.randVars = rmfield(e{1}.randVars,'genderList');
    e{1}.randVars = rmfield(e{1}.randVars,'contrastList');
    e{1}.randVars = rmfield(e{1}.randVars,'noiseList');
    %% Staircase Fix:
    % Staircase is, (NOISE1/CON2),CUE1/CUE4,PED1:3
    e{1}.stairNoise = squeeze(stimulus.staircase(1,:,:));
    e{1}.stairNoiseDual = squeeze(stimulus.dualstaircase(1,:,:));
    e{1}.stairCon = squeeze(stimulus.staircase(2,:,:));
    e{1}.stairConDual = squeeze(stimulus.dualstaircase(2,:,:));
    e{2}.stair = squeeze(stimulus.p.staircase);
    e{2}.stairDual = squeeze(stimulus.p.dualstaircase);
    %% Peripheral Re-organization
    %           f m non
    per_stim = [9 8 10];
    e{2}.randVars.tGen = zeros(size(e{2}.randVars.position));
    e{2}.randVars.tImage = zeros(size(e{2}.randVars.position));
    e{2}.correct = zeros(size(e{2}.randVars.position));
    for i = 1:length(e{2}.randVars.position)
        pos = e{2}.randVars.position(i);
        e{2}.randVars.tGen(i) = e{2}.randVars.gender{i}(pos);
        e{2}.randVars.tImage(i) = e{2}.randVars.images{i}(pos);
        if isnan(e{2}.response(i))
            e{2}.correct(i) = nan;
        else
            if e{2}.response(i) == 10 && e{2}.randVars.respond(i) == 0
                e{2}.correct(i) = 1;
            elseif e{2}.randVars.respond(i) == 1
                e{2}.correct(i) = e{2}.response(i)==per_stim(e{2}.randVars.tGen(i));
            else
                e{2}.correct(i) = 0;
            end
        end
    end
    e{2}.randVars = rmfield(e{2}.randVars,'gender');
    e{2}.randVars = rmfield(e{2}.randVars,'images');
    %% run Variables
    e{3}.runVars.blocks = stimulus.blocks.blockList;
    e{3}.runVars.dual = stimulus.dualList(end);
    e{3}.runVars.pedestals = stimulus.pedestals;
    e{3}.runVars.runNum = stimulus.counter;
    e{3}.stimulus = stimulus;
    %% Save
    save(outFile,'e');
    success = 1;
catch
    success = 0;
end