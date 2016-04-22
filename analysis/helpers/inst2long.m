function data = inst2long( inst, coh )
disp('warning: only works for coherence in stablecontrast... fix?');
%INST2LONG Convert a instances structure, returned from getInstances into a
%long form dataset.
tCount = 0; maxV = 0;
for ri = 1:length(inst)
    cinst = inst{ri}.classify.instances;
    for ii = 1:length(cinst)
        tCount = tCount + size(cinst{ii},1);
        if size(cinst{ii},2)>maxV, maxV = size(cinst{ii},2); end
    end
end
tCount = tCount*maxV; % max possible size

% setup data to the max possible size, we will remove the end later
data = zeros(tCount,5);

count = 1;
disp('Starting conversion to long form...');
disppercent(-1/length(inst));
for li = 1:length(inst)
    cinst = inst{li};

    if ~isempty(cinst.classify.instances) 
        amps = cinst.classify.instances;

%         roi_pos = find(cellfun(@(x) strcmp(x,inst.name),neural.ROIs));
        
        for ci = 1:length(amps)
            coherence = coh(ci);
            amp = amps{ci};
            for ai = 1:size(amp,1)
                for vi = 1:size(amp,2)
                    % header: ROI voxel amp repeat coherence
                    dat = [li vi amp(ai,vi) ai coherence];
                    data(count,:) = dat;
                    count = count+1;
                end
            end
        end
    end

    disppercent(li/length(inst));
end
disppercent(inf);

data = data(1:count-1,:);