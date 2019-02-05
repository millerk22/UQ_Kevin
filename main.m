%%
addpath(genpath(pwd))

load('glomo.mat')
fid = sample_fidelity(new_gt, 1, 0.1);
%% Nystrom
window = 5;
H = movmean(H,window,2);
nystrom_opt             = {};
nystrom_opt.tau         = 1/40000;
nystrom_opt.Laplacian   = 'n';
nystrom_opt.Metric      = 'Euclidean';
nystrom_opt.numsample   = 6000;
nystrom_opt.neig        = 1000;
disp('Start Nystrom')
tic
[phi, E] = nystrom(H', nystrom_opt);
clear H
%save('./data/VE_9_videos', 'nystrom_opt', 'phi', 'E');
time1_ = toc;
disp(['finishe Nystrom', num2str(time1_),'seconds']);
%%
%save('./data/HUJI/VE.mat', 'nystrom_opt', 'phi', 'E');
%%
disp('Running MBO...');
tic;
max_iter =100;
dt = 0.1;
eta = 400; 
[m1, u_hat] = mbo_multiclass(dt, eta, max_iter, phi, E, fid); 
time_ = toc;

disp(['Running MBO took ', num2str(time_),  ' s']);

%%
[~, temp] = max(m1, [], 2);
acc = sum(temp(:) == new_gt(:)) / length(new_gt);
disp(['Classification Accuracy = ', num2str(acc)]);
save('mbo_output2.mat','acc', 'temp', 'new_gt', 'nystrom_opt','eta','dt','window','fid');
