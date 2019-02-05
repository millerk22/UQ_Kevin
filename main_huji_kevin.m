
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes: Sample from the lowest 50% fidelity seems to be the best
% 181104_1 gamma = 0.05, lowest 50%, 50 iterrations to 50%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loading Data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function main_huji_kevin(suffix, niter, total_fid, top_cutoff, bottom_cutoff, cumulative_fid)

            % the first function splits the data into training and testing,
            % saving it in the /temp directory. 
ratio =  huji_sampling_60seg(1);  % ratio is the # of ground truth labels (not really ratio)
global phi E data_root
%suffix      = '_181104_1';           
data_root = './data/HUJI/';         % Location of unshuffled data
temp_root = './data/HUJI/temp/';    % Shuffled data
load([data_root,'ground_truth_huji_disney_60seg.mat']) % load the ground truth vector (in unshuffled ordering)
load([temp_root, 'train_test_index.mat'], 'train_index', 'test_index'); % these are the indices (in the unshuffled vector) that correspond to the train/test split.
fprintf("Training set: %2f%%\n", 100 * numel(train_index)/numel(ground_truth))



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nystrom 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nystrom_opt             = {};
nystrom_opt.tau         = -2;
nystrom_opt.K = 40;
nystrom_opt.Laplacian   = 'n';
nystrom_opt.Metric      = 'Euclidean';
nystrom_opt.numsample   = 400;
outstr = struct_description(nystrom_opt);
Nystrom_huji(nystrom_opt);          % Note: this Nystrom_huji function 
      % is contained in this file, finds precomputed evals/evecs if present
load([data_root,'VE_',outstr,'.mat'],'phi','E')
shuffle_index = [train_index test_index];
phi = phi(shuffle_index,:);     % phi is global variable for the eigenvectors
ground_truth = ground_truth(shuffle_index);  % Here we are shuffling the ground_truth to follow the same ordering as the shuffled data in /temp
fid_opts.type = 'random';    % We waited to load this (above) at this point to allow us to use precomputed evecs each time (independent of the shuffling done)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Down sampling?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ratio           = floor(ratio / 5);
ground_truth    = ground_truth(1:5:length(ground_truth));
phi             = phi(1:5:size(phi,1),:);



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run UQ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MBO_opts.dt = 0.05;
MBO_opts.eta = 300;
MBO_opts.neig = 400;

% define the fidelity (semi-supervision for our training)
init_fid = sample_fidelity(ground_truth, ratio/numel(ground_truth),total_fid/niter,fid_opts);
fid = init_fid;
posterior_mean  = {}; 
class_result    = {};   class_result_MBO = {};
summary_stats   = {};   summary_stats_MBO = {};
fid_perc        = [total_fid/niter];


for i = 1:niter
    fid_list = get_fid_list(fid);       % convert the MATLAB array to a list for future checks
    fprintf("Iter %d/%d -------------------------\n", i,niter);
    [posterior_mean{i}, class_result{i}]...
                            = Gibbs_huji(fid);      % perform the Gibbs sampling
    class_result_MBO{i}     = MBO_huji(fid,MBO_opts);   % perform MBO
    
    % summary stats = precision, recall, accuracy, confusion matrices
    summary_stats{i}        = get_summary_stats(fid, class_result{i}, ...
                                                ground_truth, ratio);
    summary_stats_MBO{i}    = get_summary_stats(fid, class_result_MBO{i}, ...
                                                ground_truth, ratio);
    fprintf("Mean Recall: %f\n", summary_stats_MBO{i}.recall);
    fprintf("Accuracy: %f\n", summary_stats_MBO{i}.acc);
    fprintf("Mean Precision: %f\n", summary_stats_MBO{i}.precision);
    fprintf("Fidelity percentage: %f ------------------\n", fid_perc(i));

    
%-------------------------------------------------------------------------
% Select new fidelity points!
%   We select the new fidelity points by seeing which points (in the
%   possible training set, <= ratio) have smallest variance --> selecting 
% the points in the training set that we are the least sure about..
%--------------------------------------------------------------------------

    confidence = var(posterior_mean{i},[],2);   % get the variance of the sample mean of each entry across the different classes 
    [~,idx] = sort(confidence);             % higher variance --> 
    idx = idx(idx <= ratio);                % Restrict  to the training set
    idx = idx(~ismember(idx,fid_list));     % Removing existing fidelity nodes
    fid_size = zeros(numel(fid),1);
    
    % Go through each class in the fidelity array and add more fidelity
    % points. We look at the points that aren't already fidelity nodes
    % (from exclusion above) and add a random selection of points that had 
    % low variance. We have upper and lower cutoff as percentages of the
    % total number of points in idx.
    for k = 1:numel(fid)
        idx_per_class   = idx(ground_truth(idx) == k);  % restrict to just this class in training set
        n               = sum(ground_truth(1:ratio) == k); % num in this class in training set
        num_additional_points   = ceil((total_fid/niter) * n); % number of additional points that we can add FOR THIS class
        idx_per_class           = idx_per_class(1+floor(bottom_cutoff * length(idx_per_class)) ... % 
                                                :ceil(top_cutoff * length(idx_per_class))); 
        idx_per_class           = idx_per_class(randperm(length(idx_per_class)));
        % shuffle;
        fid{k}                  = union(fid{k}, ...
                                        idx_per_class(1:num_additional_points));
        fid_size(k) = numel(fid{k});    % change the fidelity size entries to reflect number of points now in fidelity for this class
    end
    fid_perc(i+1) = sum(fid_size)/ratio;        % update fidelity percentage (% of training set that we are currently looking at)
    
    % Save the results of this iteration
    save(['iter_UQ_results_huji', suffix, '.mat'],'posterior_mean', 'class_result', ...
    'fid_perc', 'summary_stats', 'class_result_MBO', 'summary_stats_MBO')
end



fid = init_fid;         % Reset the fidelity to be the original fidelity points 
                    % (the part of the training data chosen as initial
                    % fidelity)
class_result_noniter = {}; posterior_mean_noniter = {};
summary_stats_noniter = {}; summary_stats_noniter_MBO = {};
class_result_noniter_MBO = {};
fprintf("Start non-iterative UQ --------------------\n")
% Non-iterative UQ..  We don't iteratively add the least sure training data
% in each time. We either:
%     1 - sample the additional points of the training points, adding them
%     to the fidelity at next iteration. This is like before, BUT not
%     adding based on the confidence score.
%     2 - sample a new portion of the fidelity points (cumulative_fid =
%     false)


for i=1:niter
    fid_list = get_fid_list(fid);
    [posterior_mean_noniter{i}, class_result_noniter{i}] = Gibbs_huji(fid);
    class_result_noniter_MBO{i} = MBO_huji(fid, MBO_opts);
    summary_stats_noniter{i} = get_summary_stats(fid, ...
        class_result_noniter{i}, ground_truth, ratio);
    summary_stats_noniter_MBO{i} = get_summary_stats(fid, ...
        class_result_noniter_MBO{i}, ground_truth, ratio);
    fprintf("Mean Recall: %f\n", summary_stats_noniter_MBO{i}.recall);
    fprintf("Accuracy: %f\n", summary_stats_noniter_MBO{i}.acc);
    fprintf("Mean Precision: %f\n", summary_stats_noniter_MBO{i}.precision);
    fprintf("Fidelity percentage: %f ------------------\n", fid_perc(i));
    fid_size = zeros(numel(fid),1);
    if cumulative_fid
        idx = 1:ratio;
        idx = idx(~ismember(idx,fid_list));     % Removing existing  fidelity nodes
        fid_size = zeros(numel(fid),1);
        for k = 1:numel(fid)
            idx_per_class   = idx(ground_truth(idx) == k);
            n               = sum(ground_truth(1:ratio) == k);
            num_additional_points   = ceil((total_fid/niter) * n);
            idx_per_class = idx_per_class(randperm(length(idx_per_class)));
        % shuffle
            fid{k} = union(fid{k}, ...
                   idx_per_class(1:num_additional_points));
            fid_size(k) = numel(fid{k});
        end
    else
        fid  =  sample_fidelity(ground_truth, ...
            ratio/numel(ground_truth),fid_perc(i+1),fid_opts);
        for k = 1:numel(fid)
            fid_size(k) = numel(fid{k});
        end
    end
    fid_perc(i+1) = sum(fid_size)/ratio;
end
save(['noniter_UQ_results_huji', suffix, '.mat'],'posterior_mean_noniter',...
    'class_result_noniter', 'summary_stats_noniter', ...
    'class_result_noniter_MBO','summary_stats_noniter_MBO', 'fid_perc')
end

function [temp] = MBO_huji(fid,MBO_opt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MBO_opt:
%           dt      step size
%           eta     weights before fidelity
%           neig    number of eigenvalues/vectors to use
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global phi E
    max_iter = 100; 
    disp('Running MBO...');
    tic;
    [m1, ~] = mbo_multiclass(MBO_opt.dt, MBO_opt.eta, max_iter, ...
                            phi(:,1:MBO_opt.neig), E(1:MBO_opt.neig), fid); 
    time_ = toc;
    disp(['Running MBO took ', num2str(time_),  ' s']);
    [~, temp] = max(m1, [], 2);
end


function [m1, sol] = Gibbs_huji(fid)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MBO_opt:
%           dt      step size
%           eta     weights before fidelity
%           neig    number of eigenvalues/vectors to use
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global phi E
    max_iter = 1.5e4; 
    gamma = 0.05;
    disp('Running Gibbs...');
    tic;
    [m1, ~] = mcmc_gibbs_probit_multiclass(max_iter, gamma,  ...
                            phi, E,fid); 
    time_ = toc;
    [~,sol] = max(m1,[],2);
    disp(['Running Gibbs took ', num2str(time_),  ' s']);
end


function Nystrom_huji(nystrom_opt)
global data_root
    outstr = struct_description(nystrom_opt);
    fp = [data_root, 'VE_',outstr,'.mat'];
    if exist(fp,'file')
       fprintf('Found precomputed eigenvectors/values\n');
       return
    end
    load([data_root,'H_huji_disney_60seg.mat'],'H');
    H = movmean(H,20,2);  % why do we do a moving mean here? To reduce noise?
    disp("Performing Nystrom extension")
    tic;
    [phi, E] = nystrom(H', nystrom_opt);
    t = toc;
    fprintf("Nystrom takes %f seconds\n", t);
    clear H 
    save(fp, 'nystrom_opt', 'phi', 'E');
end
