% We now have, in /f and /m, two datasets of faces. Our next job is to sort
% these to determine whether any are "ambiguous" in gender. 

cFolder = '~/proj/att_awe/images/brazil_faces/';

gen = {'m/','f/'};

%% Load face data
data.images = zeros(193,162,200);
data.fname = {};
data.gender = [];
mcount = 0;
fcount = 0;


counter = 1;
for gi = 1:length(gen)
    gender = gen{gi};
    loc = [cFolder gender];
    imgFiles = dir(strcat(loc,'*jpg'));
    for imgNum = 1:length(imgFiles)
        im = imread([loc imgFiles(imgNum).name]);
        data.images(:,:,counter) = im;
        data.fname{counter} = imgFiles(imgNum).name;
        data.gender(counter) = gi;
        counter = counter + 1;
    end
end

%% Load data holder
cdata = struct;
cdata.name = 'dan';
cdata.gender = [];
cdata.keep = ones(1,200);

% Important attributes:
%   Estimated Gender
%   Attractiveness
%   Eyes -> Camera

%%
figure(1)
colormap('gray')

%% Gender
randlist = randperm(200);
sidx = 1;
for i = sidx:200
    ci = randlist(i);
    imagesc(data.images(:,:,ci));
    cdata.gender(ci) = str2num(input('Male, Not 100% Sure, Female? ','s'));
end

%% Check for accidents
list = 1:200;
mlist = list(cdata.gender==1);
mlist = mlist(mlist>100);
disp('Listed as MALE');
for i = 1:length(mlist)
    ci = mlist(i);
    imagesc(data.images(:,:,ci));
    cdata.keep(ci) = str2num(input('Keep?: ','s'));
end
flist = list(cdata.gender==3);
flist = flist(flist<100);
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

%% Remove files
for gi = 1:length(gen)
    gender = gen{gi};
    loc = [cFolder gender];
    for i = 1:length(cdata.keep)
        if ~cdata.keep(i)
            delete([loc data.fname{i}]);
        end
    end
end
