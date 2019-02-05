function [ y_mean,iter_stats ] = mcmc_gibbs_probit_multiclass_kevin(max_iter, gamma, V, E, fid, varargin)
%--------------------------------------------------------------------------
% Author: Xiyang Luo <xylmath@gmail.com> , UCLA
% ******Comments for the sake of Kevin Miller... I believe the description of
% Inputs etc here is out of date.. (first input is not beta --> beta
% corresponds to the pCN algorithm). 


% This file is part of the diffuse-interface graph algorithm code.
% There are currently no licenses.
%
%--------------------------------------------------------------------------
%  Description :
%       Gibbs Sampler
%  Model :
%           P(u) = exp(- J(u)/T), where J(u) = .5*<u, Lu> + Phi(u)
%  Proposal :
%           (1) x' = sqrt(1-beta^2) * x + beta L^-1/2 w
%           (2) acc = min{1, likelihood(x') / likelihood(x)}
%
%  Input :
%       beta: scalar, proposal variance for MCMC
%       max_iter: maximum of iteration for MCMC
%       Phi_eval: likelihood function handle
%       V: eigenvectors of L. (used to compute iteration stats only).
%       E: normalized eigenvalues.
%       fid : NClassx1 cell array, indices of fidelity points for each class
%       varargin: opt: options (See Description in
%                               make_pcn_default_options.m)
%
%  Output :
%       y_mean : Expectation of the label prediction y = S(u)
%       iter_stats : stats for the MCMC iterations. Details below:
%                           vvv
%       y_mean_rec : record of the mean of the y-prediction.
%       u_rec : record of the current sample of u.
%       u_hat_rec : record of u_hat, i.e., coefficient on eigenvectors
%       mcmc_cum_prob : cummulative acceptance probability of the MCMC.
%       acceptance_prob : raw acceptance probability for each iteration
% -------------------------------------------------------------------------



% default options
if(size(varargin) == 0)
    opt = {};
else
    opt = varargin{1};
end
opt = make_pcn_default_options_binary(opt);  % binary and multiclass default options are exactly
            % the same
%opt = make_pcn_default_options_multiclass(opt);


% initialize variables
if opt.BurnIn > 0
    rec_length = max_iter - opt.BurnIn;
else
    rec_length = max_iter;
end
N = size(V,1);
Neig = size(V,2);
NClass = numel(fid);
trandn_obj = trandn_multiclass(NClass); % this object will be used to 
        % generate samples from a truncated normal sampler for multiclass


Nf = zeros(NClass, 1); % Nf is number of fidelity, for each class
for i = 1 : NClass
    Nf(i) = numel(fid{i}); 
end
num_fid = sum(Nf); 
%u_hat_rec=zeros(Neig, NClass, rec_length);  % from previous copy?
y_mean=zeros(N,NClass);  % mean?
u_mean = y_mean;
%v_rec=zeros(num_fid, NClass, rec_length);  % from previous copy?


if opt.isrec_u  % if want to record results, then we instatiate variables for that
    u_rec=zeros(N, NClass, rec_length);
    y_mean_rec=zeros(N, NClass, rec_length);
    u_mean_rec=zeros(N, NClass, rec_length);
end

% sample in batch a Gaussian Noise.
z_hat_all=randn(Neig, NClass, max_iter);


% initialize u
u=zeros(N, NClass);
if(ischar(opt.Init))
    if strcmp(opt.Init, 'fid') % if have fidelity, set to corr one-hot vectors
        for cl = 1:NClass
            for i = 1:numel(fid{cl})
                u(fid{cl}(i), cl) = 1;
            end
        end
    end
    if strcmp(opt.Init, 'rand') % random initialization
        u = multiclass_threshold(randn(N, NClass));
    end
else
    u = u.Init;     % I believe this is wrong.. this probably should be opt.Init, the option
end                 % in the case that opt.Init is an actual initial distribution?..




% permute indices for better vectorization -- put fidelity at beginning
node_index = zeros(N, 1);
count = 1; 
for i = 1 : NClass
    node_index(count : Nf(i) + count - 1) = fid{i}; 
    count = count + Nf(i); 
end
remain_index = 1:N; 
for i = 1 : NClass
    remain_index = setdiff(remain_index, fid{i});
end
node_index(count : end) = remain_index(:);
[~, reverse_node_index] = sort(node_index, 'ascend');
V = V(node_index, :);
u = u(node_index, :);



% Precompute E^-1/2 and force samples to lie in the orthogonal complement
% of e^1 = (1, 1, ... 1)
gamma_sq = gamma * gamma;

% This is a preparatory calculation for the Gibbs sampling of the fidelity
% terms. See the multi-class paper, p. 3 left column. We have a
% transformation from u to xi, V_KJ here is the calculation HQ' in the
% paper. Then we are forming 

V_KJ = V (1 : num_fid, :); 
A_KJ = V_KJ' * V_KJ;       % the fidelity indices' corresponding evecs outer product 
P_KJ = diag(E) + 1/gamma_sq * A_KJ;  
P_KJ = .5 * (P_KJ + P_KJ');  % symmetrize, isn't P_KJ already symmetric though?
[Q_KJ, S_KJ] = eig(P_KJ);
S_KJ = diag(S_KJ);
invsq_skj = 1./sqrt(S_KJ); 


% Main iteration.
for k=1:max_iter
    ind = k;
  
    % Sample v ~ P(v|u)
    zk_hat = z_hat_all(:, :, ind);      %rv's for the KL expansion
    v = zeros(num_fid, NClass);
    count = 1; 
    
    % get truncated normal samples for Gibbs sampling updates 
    for i = 1 : NClass
        v(count : count + Nf(i) - 1, :) = trandn_obj.gen_samples(u(count : count + Nf(i) - 1, :), gamma, i); 
        count = count + Nf(i); 
    end
    
    
    % Sample u ~ P(u|v) via the KL expansion u = Phi xi
    temp = V_KJ'*v/gamma_sq;
    temp = Q_KJ'*temp;
    temp = bsxfun(@rdivide, temp, S_KJ); % P^-1 v
    % temp = temp./S_KJ;
    m_hat = Q_KJ * temp;
    u_hat = Q_KJ * (bsxfun(@times, invsq_skj, zk_hat)) + m_hat;
    % u_hat = Q_KJ * (invsq_skj .* zk_hat) + m_hat;
    u = V * u_hat; 

    
    % updating and recording all variables
    k_rec = k - opt.BurnIn;
    if k_rec > 0
        if k_rec == 1
            y_mean=multiclass_threshold(u);
            u_mean = u;
        else
            y_mean=((k_rec-1)*y_mean+multiclass_threshold(u))/k_rec;
            u_mean = ((k_rec-1)*u_mean + u)/k_rec;
        end
        if opt.isrec_u
            u_rec(:, :, k_rec)=u;
            y_mean_rec(:, :, k_rec)=y_mean;
            u_mean_rec(:, :, k_rec) = u_mean;
        end
        %u_hat_rec(:, :, k_rec)=u_hat;
        %v_rec(:, :, k_rec) = v; 
        % Display Verbose
        if opt.Verbose && mod(k_rec, opt.VerboseWindowLength) == 0
            disp(['Iteration ', num2str(k)]);
        end
    end
end

% Reverse the permutation of everything
y_mean = y_mean(reverse_node_index, :, :);
u_mean = u_mean(reverse_node_index,:,:);

% record iteration results to iter_stats.
if opt.isrec_u
    u_rec = u_rec(reverse_node_index, :, :);
    y_mean_rec = y_mean_rec(reverse_node_index, :);
    iter_stats.u_rec = u_rec;
    iter_stats.y_mean_rec = y_mean_rec;
end
%iter_stats.u_hat_rec = u_hat_rec;
%iter_stats.v_rec = v_rec; 
iter_stats.u_mean = u_mean;
end




