function fullData = loadFullData( )


if isfile(fFile)
    disp('(fullData) Loading...');
    load(fFile);
else
    disp('(fullData) Generating New...');
    fullData = struct;
end