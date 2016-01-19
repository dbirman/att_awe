function cleaned_files = cc_clean( files )
%CC_CLEAN Remove voxels with no response or ridiculous responses
%
% Restrict all voxels to:
% < 20 amplitude
% > -20 amplitude

cleaned_files = cell(1,length(files));

for fi = 1:length(files)
    fname = files{fi};
    load(fname)
    % ldata and rdata get loaded
    % remove spikes
    ldata = ldata(~any(ldata<-100,2),:);
    ldata = ldata(~any(ldata>100,2),:);
    % remove more crap
    lkeep = logical(logical(~any(ldata<-100,1)) .* logical(~any(ldata>100,1)) .* ~any(isnan(ldata),1));
    ldata = ldata(:,lkeep);
    
    % remove spikes
    rdata = rdata(~any(rdata<-100,2),:);
    rdata = rdata(~any(rdata>100,2),:);
    % remove more crap
    rkeep = logical(logical(~any(rdata<-100,1)) .* logical(~any(rdata>100,1)) .* ~any(isnan(rdata),1));
    rdata = rdata(:,rkeep);
    
    % z-score
    hold = ldata(:,5:end);
    ldata(:,5:end) = (hold-repmat(mean(hold,1),size(hold,1),1))./repmat(std(hold,1),size(hold,1),1);
    
    fout = strrep(fname,'dataset','cleanset');
    cleaned_files{fi} = fout;
    save(fout,'ldata','rdata');
end

end

