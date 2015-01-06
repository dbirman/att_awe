% We now have, in /f and /m, two datasets of faces. Our next job is to sort
% these to determine whether any are "ambiguous" in gender. 

cFolder = '~/proj/att_awe/images/real_faces2/';

gen = {'m/','f/'};

%% Load face data
data.images = zeros(100,100,182);
data.fname = {};
data.gender = [];
mcount = 0;
fcount = 0;

mask_size = [100 100];
mask = repmat([zeros(1,75),ones(1,100),zeros(1,75)],100,1);
mask = [zeros(90,250);mask;zeros(60,250)];

counter = 1;
for gi = 1:length(gen)
    gender = gen{gi};
    loc = [cFolder gender];
    imgFiles = dir(strcat(loc,'*jpg'));
    for imgNum = 1:length(imgFiles)
        im = imread([loc imgFiles(imgNum).name]);
        im = reshape(im(mask==1),mask_size(1),mask_size(2));
        data.images(:,:,counter) = im;
        data.fname{counter} = imgFiles(imgNum).name;
        data.gender(counter) = gi;
        counter = counter + 1;
    end
end

%% Load data holder
cdata = struct;
cdata.name = 'dan';
cdata.number = 1;
cdata.gender = [];
cdata.attract = [];
cdata.age = [];
cdata.keep = ones(1,182);

% Important attributes:
%   Estimated Gender
%   Attractiveness
%   Eyes -> Camera

%%
figure(1)
colormap('gray')

%% Gender
randlist = randperm(182);
sidx = 1;
for i = sidx:182
    ci = randlist(i);
    imagesc(data.images(:,:,ci));
    cdata.gender(ci) = str2num(input('Male, Not 100% Sure, Female? ','s'));
end

%% Check for accidents
list = 1:182;
mlist = list(cdata.gender==1);
mlist = mlist(mlist>93);
disp('Listed as MALE');
for i = 1:length(mlist)
    ci = mlist(i);
    imagesc(data.images(:,:,ci));
    cdata.keep(ci) = str2num(input('Keep?: ','s'));
end
flist = list(cdata.gender==3);
flist = flist(flist<94);
disp('Listed as FEMALE');
for i = 1:length(flist)
    ci = flist(i);
    imagesc(data.images(:,:,ci));
    cdata.keep(ci) = str2num(input('Keep?: ','s'));
end
alist = list(cdata.gender==2);
disp('Listed as AMBIGUOUS');
for i = 1:length(alist)
    ci = alist(i);
    imagesc(data.images(:,:,ci));
    cdata.keep(ci) = str2num(input('Keep?: ','s'));
end

%% Remove the no-keep
randlist = 1:182;
randlist = randlist(cdata.keep==1);
randlist = randlist(randperm(length(randlist)));

%% Remove files
for i = 1:length(cdata.keep)
    if ~cdata.keep(i)
        delete([loc data.fname{i}]);
    end
end

%% Attractive
sidx = 1;
for i = sidx:182
    ci = randlist(i);
    imagesc(data.images(:,:,ci));
    cdata.attract(ci) = str2num(input('Attractive? 1 = Not Very, 5 = Very: ','s'));
end

%% Age
randlist = randperm(182);
sidx = 1;
for i = sidx:182
    ci = randlist(i);
    imagesc(data.images(:,:,ci));
    cdata.age(ci) = str2num(input('Age? Input Estimate in 5s: ','s'));
end

%% Find a set of 80 male/female files with an even mean age of 40\

cdata.mattract = cdata.attract(1:93);
cdata.fattract = cdata.attract(94:end);

mnums = logical([ones(1,75) zeros(1,18)]); % we will drop 8 male files
fnums = logical([ones(1,75) zeros(1,14)]); % we will drop 4 female files

%%
minN = 50; best = [];
while(mean(cdata.mattract(mnums)) > 42 || mean(cdata.mattract(mnums)) < 38)
    if mean(cdata.mattract(mnums)) < minN
        minN = mean(cdata.mattract(mnums));
        best = mnums;
        disp(minN);
    end
    mnums = mnums(randperm(length(mnums)));
end
bestM = best;
%%
maxN = 30; best = [];
while(mean(cdata.fattract(fnums)) > minN+1 || mean(cdata.fattract(fnums)) < minN-1)
    if mean(cdata.fattract(fnums)) > maxN 
        maxN = mean(cdata.fattract(fnums));
        best = fnums;
        disp(maxN);
    end
    fnums = fnums(randperm(length(fnums)));
end
bestF = best;

%% Copy all the good files

cFolderM = '~/proj/att_awe/images/real_faces2/m/';
oFolderM = '~/proj/att_awe/images/real_faces2_agem/m/';
cFolderF = '~/proj/att_awe/images/real_faces2/f/';
oFolderF = '~/proj/att_awe/images/real_faces2_agem/f/';

for i = 1:93
    if bestM(i)==1
        copyfile([cFolderM data.fname{i}],[oFolderM data.fname{i}]);
    end
end

for i = 94:182
    if bestF(i-93)==1
        copyfile([cFolderF data.fname{i}],[oFolderF data.fname{i}]);
    end
end