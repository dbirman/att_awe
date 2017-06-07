%% Simulate poisson with correlations

% Simulate 36 poisson processes representing gaussian tuning functions. For
% each process assume the underlying population consists of 50 units


popmean = 0;
popstd = 25;
popn = 10;

structure = randi(5,1,popn);

angles = 0:10:350;
contrasts = 10:10:100;

data = zeros(length(contrasts),36,36,popn);

for ci = 1:length(contrasts)
    con = contrasts(ci);
    for pi = 1:length(angles)
        population = angles(pi);
        for ti = 1:length(angles)
            theta = angles(ti);
            % this population prefers angle population, find my firing rate
            spkrate = 10+con*normpdf(theta-population,0,popstd);
            data(ci,pi,ti,:) = structure .* poissrnd(spkrate,1,popn);
        end
    end
end

%% Compute the noise across any particular tuning
figure; hold on
for ci = 1:length(contrasts)
    con = contrasts(ci);
    cdat = squeeze(data(ci,:,:,:));
    
    mu = mean(cdat(:));
    sd = std(cdat(:));
    
    plot(con,mu,'*r');
    plot(con,sd,'*b');
    
end
legend({'Mean','SD'});