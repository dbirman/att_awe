function cc_slope(type, ocon, ocoh )
%%
% get slope over ROIs per subject
% slope = zeros(size(ocon,1),1);

% ACROSS rois
% for ni = 1:size(ocon,1)
%     y = ocoh(ni,:)';
%     x = ocon(ni,:)';
%     x = [ones(size(x)) x];
%     b = x\y;
%     slope(ni) = b(2);
% end

% ACROSS SUBJECTs (definitely wrong)
for ri = 1:size(ocon,2)
    y = ocoh(:,ri);
    x = ocon(:,ri);
    x = [ones(size(x)) x];
    b = x\y;
    slope(ri) = b(2);
end

% DIFFERENCE
oc = ocon(:);
om = ocoh(:);
if any(oc>100)
    om = om(oc<100);
    oc = oc(oc<100);
end
if any(om>100)
    oc = oc(om<100);
    om = om(om<100);
end
if any(oc<-100)
    om = om(oc>-100);
    oc = oc(oc>-100);
end
if any(om<-100)
    oc = oc(om>-100);
    om = om(om>-100);
end
    
ci = bootci(10000,@mean,oc);
ci
mean(ci)

ci = bootci(10000,@mean,om);
ci
mean(ci)

slope = oc(:)-om(:);


ci = bootci(10000,@mean,slope);

disp(sprintf('%s slope = %02.2f 95%% CI [%02.2f %02.2f]',type,mean(ci),ci(1),ci(2)));
