
addpath(genpath(pwd))
load('glomo.mat','new_gt','H');
batch_size = 30000;
temp = {};
%new_gt = {};
%ugt = {};
%total_mask = zeros(floor(length(gt)/batch_size),1);
for batch_idx = 1:floor(length(new_gt)/batch_size)
    idx = ((batch_idx -1)* batch_size + 1): (batch_idx * batch_size);
    data = H(:, idx);
    %data = rand(size(data));
    gt_data = new_gt(idx);
    fid = sample_fidelity(gt_data, 1,0.1);
    smooth_data = movmean(data,5,2);
    %smooth_data = data;
    tic
    nystrom_opt             = {};
    nystrom_opt.tau         = -2;
    nystrom_opt.K = 100;
    nystrom_opt.Laplacian   = 'n';
    nystrom_opt.Metric      = 'Euclidean';
    nystrom_opt.numsample   = 3000;
    nystrom_opt.neig        = 100;
    disp('Start Nystrom')
    tic
    [phi, E] = nystrom(smooth_data', nystrom_opt);
    %save('./data/VE_9_videos', 'nystrom_opt', 'phi', 'E');
    time1_ = toc;
    disp(['finishe Nystrom', num2str(time1_),'seconds']);
    dt = 0.1;
    eta =400.0;
    max_iter = 100;
    %%
    [u, u_hat] = mbo_multiclass(dt, eta, max_iter, phi, E, fid);
    [~,temp{batch_idx}] = max(u,[],2);
    temp{batch_idx} = recover_batch_label(gt_data, temp{batch_idx});
    [ C ] = confusion_matrix(gt_data, Nind2vec(temp{batch_idx} ));
    C_total{batch_idx} = C(unique(gt_data),unique(gt_data));
    save('./results/mbatch_glomo_neig100_nsample_3000_3K.mat', 'new_gt', 'temp', 'dt', 'eta', 'nystrom_opt','C_total');

    %%
end

save('./results/mbatch_glomo_neig100_nsample_3000_3K.mat', 'new_gt', 'temp', 'dt', 'eta', 'nystrom_opt','C_total');
%mean(acc)
