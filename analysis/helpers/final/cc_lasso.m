function outputs = cc_lasso( neural, sid, datasets )
%CC_LASSO Summary of this function goes here
%   Detailed explanation goes here

outputs = {};

for fi = 1:length(datasets)
    disp('************************************');
    disp(sprintf('Computing CV Lasso for Folder %s',neural.folders{fi}));
    disp('************************************');
    load(datasets{fi});
    % loaded ldata rdata
    out = struct;
    out.ldata = ldata;
    out.rdata = rdata;
    
    out.left = struct;
    out.right = struct;
    
    %% LEFT
    out.left = cc_lasso_helper(ldata);
    out.right = cc_lasso_helper(rdata);
    
    %% Figure
    figfol = '~/proj/att_awe/analysis/figures/decode/';
    
    tasks = {'coh','con'};
    for ti = 1:2
        task = tasks{ti};
        figfile = fullfile(figfol,sprintf('%s_%s_attending_%s',sid,neural.folders{fi},task));
        figure, hold on
        title(sprintf('Attending %s',task));
        dat = out.left.(task);
        for foldi = 1:10
            plot(dat.true_con{foldi},dat.pred_con{foldi},'r*');
            plot(dat.true_coh{foldi},dat.pred_coh{foldi},'b*');
        end
        dat = out.right.(task);
        for foldi = 1:10
            plot(dat.true_con{foldi},dat.pred_con{foldi},'r*');
            plot(dat.true_coh{foldi},dat.pred_coh{foldi},'b*');
        end
        legend({'Decoding Contrast','Decoding Coherence'});
        xlabel('True Value');
        ylabel('Predicted Value');
        print(figfile,'-dpdf');
    end
    
    outputs{end+1} = out;

end

