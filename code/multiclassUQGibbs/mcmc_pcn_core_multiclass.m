function [ y_mean1,iter_stats ] = mcmc_pcn_core_multiclass(beta, max_iter, Phi_eval, V, E, fid, varargin)
%--------------------------------------------------------------------------
% Author: Xiyang Luo <xylmath@gmail.com> , UCLA
%
% This file is part of the diffuse-interface graph algorithm code.
% There are currently no licenses.
%
%--------------------------------------------------------------------------
%  Description :
%       PCN algorithm main framework for multiclass classification.
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
%       other quantities of interest. 
% -------------------------------------------------------------------------

% default options
if(size(varargin) == 0)
    opt = {};
else
    opt = varargin{1};
end
opt = make_pcn_default_options_multiclass(opt);
NumAdditionalFunc = numel(opt.AdditionalAvgFunc); 
if NumAdditionalFunc >0
    iter_stats.additional_avg_func = cell(NumAdditionalFunc, 1); 
end

% initialize variables
if opt.BurnIn > 0
    rec_length = max_iter - opt.BurnIn;
else
    rec_length = max_iter;
end
N = size(V,1);
Neig = size(E,1);
NClass = numel(fid);
acceptance_prob=zeros(rec_length,NClass);
mcmc_cum_prob=zeros(rec_length,NClass);
u_hat_mean_rec=zeros(Neig, NClass, rec_length);
u_hat_rec=zeros(Neig, NClass, rec_length);
y_mean1=zeros(N, NClass);
u_hat_mean=zeros(Neig, NClass);
if opt.isrec_u
    u_rec=zeros(N, NClass, rec_length);
    y_mean_rec=zeros(N, NClass, rec_length);
end
yl = 1:NClass;

%sample in batch a Gaussian Noise.
if opt.RandomBatchSize == -1
    z_hat_all=randn(Neig, NClass, max_iter);
    if opt.LambdaBar < Inf %% -- not used
        z_rand_all = randn(N, NClass, max_iter);
        z_perp_all = zeros(N, NClass, max_iter);
        for j = 1:NClass
            z_perp_all(:, j, :) = squeeze(z_rand_all(:, j, :)) - V * (V'*squeeze(z_rand_all(:, j, :))) ;
        end
        clear z_rand_all;
    end
else %% -- not used
    Bchsz = opt.RandomBatchSize;
end

%% initialized
u=zeros(N, NClass);
if(ischar(opt.Init))
    if strcmp(opt.Init, 'fid')
        for cl = 1:NClass
            for i = 1:numel(fid{cl})
                u(fid{cl}(i), cl) = 1;
            end
        end
    end
    if strcmp(opt.Init, 'rand')
        u = multiclass_threshold(randn(N, NClass));
    end
else
    u = u.Init;
end
%%
%Precompute E^-1/2 and force samples to lie in the orthogonal complement
%of e^1 = (1, 1, ... 1)
if opt.ZeroMean
    E(1) = 1;
    inv_sqrt_E = 1./sqrt(E);
    inv_sqrt_E(1) = 0;
    E(1) = 0;
else
    inv_sqrt_E = 1./sqrt(E);
end
%% -- not used
if ischar(opt.LambdaBar) && strcmp(opt.LambdaBar, 'auto')
    opt.LambdaBar = max(E(2:end));
end
if opt.LambdaBar < Inf  % # TODO: Add in auto option here
    inv_sqrt_lam_bar = 1.0/sqrt(opt.LambdaBar);
end
%%
%normalize constant
c = sum(inv_sqrt_E.^2);
if opt.LambdaBar < Inf  
    c = c + opt.LambdaBar*(N-Neig);
end
sqrt_c = sqrt(N*1/c);

Neig = length(E);

%% Main iteration.


for k=1:max_iter
    k
    for m = 1:NClass %modified 
        zk_hat = z_hat_all(:,m,k); 
        v = u;
        v(:,m) = sqrt(1-beta^2)*u(:,m) + sqrt_c*beta*V*(zk_hat.*inv_sqrt_E);
        logpu = Phi_eval(u, fid, yl);
        logpv = Phi_eval(v, fid, yl);
        acc = min(1, exp(logpu - logpv));
        % accep-reject step
        r=rand;
        if acc>r
            u = v;
        end        
        %record the acc_prob for different class
        k_rec = k - opt.BurnIn;
        if k_rec >0
            acceptance_prob(k_rec, m) = acc;
            if m==NClass
                if k_rec == 1
                    mcmc_cum_prob(k_rec,:) = acceptance_prob(k_rec,:);
                    y_mean1=multiclass_threshold(u);
                    y_mean0=1;
                else
                    mcmc_cum_prob(k_rec,:)=((k_rec - 1)*mcmc_cum_prob(k_rec - 1,:)+acceptance_prob(k_rec,:))/(k_rec);
                    y_mean1=((k_rec-1)*y_mean1+multiclass_threshold(u))/k_rec;
                    if NumAdditionalFunc > 0
                       for l = 1 : NumAdditionalFunc
                           avf = opt.AdditionalAvgFunc{l}; 
                           iter_stats.additional_avg_func{l} = ((k_rec - 1)*iter_stats.additional_avg_func{l} + avf(u))/(k_rec); 
                       end
                    end
                end
                if opt.isrec_u
                    u_rec(:, :, k_rec)=u;
                    y_mean_rec(:, :, k_rec)=y_mean;
                end
                u_hat_rec(:, :, k_rec)=V' * u;
                u_hat_mean=((k_rec-1)*u_hat_mean+V' * u)/k_rec;
                u_hat_mean_rec(:, :, k_rec)=u_hat_mean;
                % Display Verbose
                if opt.Verbose && mod(k_rec, opt.VerboseWindowLength) == 0
                    disp(['Iteration ', num2str(k), ': ', 'Cumulative Acc Prob = ', num2str(mcmc_cum_prob(k_rec,:))]);
                end
            end
        end
    end
    

    %terminate when y_mean converge
    if (k_rec>4000 && (sum(sum((y_mean0-y_mean1).^2))/sum(sum(y_mean0.^2)))<1e-10)
       break;
    else
        y_mean0 = y_mean1;
    end
       
end

% record iteration results to iter_stats.
if opt.isrec_u
    iter_stats.u_rec = u_rec;
    iter_stats.y_mean_rec = y_mean_rec;
end
iter_stats.mcmc_cum_prob = mcmc_cum_prob;
iter_stats.acceptance_prob = acceptance_prob;
iter_stats.eigs_cum_acc_prob = eigs_cum_acc_prob;
iter_stats.Num_eigs = Num_eigs;
iter_stats.k = k_rec;
end




