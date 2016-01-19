function out = cc_lasso_helper( data )
%CC_LASSO_HELPER Summary of this function goes here
%   Detailed explanation goes here
out = struct;

obs = size(data,1);

folds = 10;

tasks = {'coh','con'};

idxs = floor([0 linspace(1,obs,folds+1)]);



disp('Computing Folds--this could take a while');
for task = 1:2
    taskname = sprintf('%s',tasks{task});
    
    atrue_con = cell(1,folds);
    atrue_coh = cell(1,folds);
    apred_con = cell(1,folds);
    apred_coh = cell(1,folds);
    acon_loss = zeros(1,folds);
    acoh_loss = zeros(1,folds); 
    parfor i = 2:folds+1
        test_idxs = idxs(i):idxs(i+1);
        train_idxs = [1:idxs(i)-1 idxs(i+1)+1:obs];
        
        test = data(test_idxs,:);
        train = data(train_idxs,:);
        % con - coh - task - side - vox1 vox2 ... voxN
        
        % select subset of data for this task
        ctrain = select(train,3,task);
        ctest = select(test,3,task);
        X = ctrain(:,5:end);
        % train model
        [B_con, stats_con] = lasso(X,ctrain(:,1),'CV',10);
        [B_coh, stats_coh] = lasso(X,ctrain(:,2),'CV',10);
        % test model
        true_con = ctest(:,1);
        true_coh = ctest(:,2);
        test_X = ctest(:,5:end);
        pred_con = test_X*B_con(:,stats_con.IndexMinMSE);
        pred_coh = test_X*B_coh(:,stats_coh.IndexMinMSE);
        
        con_loss = sum((true_con-pred_con).^2);
        coh_loss = sum((true_coh-pred_coh).^2);
        
        atrue_con{i-1} = true_con;
        atrue_coh{i-1} = true_coh;
        apred_con{i-1} = pred_con;
        apred_coh{i-1} = pred_coh;
        
        acon_loss(i-1) = con_loss;
        acoh_loss(i-1) = coh_loss;
        disp('A fold finished (10 total)');
    end
    out.(taskname).true_con = atrue_con;
    out.(taskname).true_coh = atrue_coh;
    out.(taskname).pred_con = apred_con;
    out.(taskname).pred_coh = apred_coh;
    out.(taskname).con_loss = acon_loss;
    out.(taskname).coh_loss = acoh_loss;
end

function dat = select(dat,col,val)

keep = dat(:,col)==val;
dat = dat(keep,:);