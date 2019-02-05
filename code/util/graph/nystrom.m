function [phi, E] = nystrom(data, opt)

%--------------------------------------------------------------------------
% Author: Xiyang Luo <xylmath@gmail.com> , UCLA 
%
% This file is part of the diffuse-interface graph algorithm code. 
% There are currently no licenses. 
%
%--------------------------------------------------------------------------
%  Description : Nystrom Extension for Normalized Laplacian
%
%  Input : 
%       data = (x1T;x2T;...xmT), where individual points are rows
%       opt. tau : kernel width for metric
%       opt. numsample : number of samples
%       opt. Metric : 'Euclidean' or 'Cosine'
%       opt. Laplacian = 'n' or 'u'
%       opt. fid    : make sure to include these nodes as fidelity
%
%  Output : 
%       phi : (phi_1,...phi_n), coloumn eigenvectors
%       E   : Array of eigenvalues (increasing)
% -------------------------------------------------------------------------

opt = check_opt(opt); 
tau = opt.tau;
num_samples = opt.numsample;
neig = opt.numsample;

% randomly select samples
num_rows = size(data, 1);
if ~isfield(opt, 'fid')
    permed_index = randperm(num_rows);
    sampled = permed_index(1:num_samples);
    sample_data = data(sampled, :);
    other = setdiff(1:num_rows, sampled);
    other_data = data(other,:);
    permed_index = [sampled, other];
else
    fid = opt.fid;
    if (numel(fid) > num_samples)
        sample_data = data(fid,:);
        other = setdiff(1:num_rows, fid);
        other_data = data(other,:);
        permed_index = [fid, other];
        num_samples = numel(fid);
    else
        permed_index = randperm(num_rows);
        sampled = permed_index(1:num_samples);
        sampled = union(sampled, fid);
        sample_data = data(sampled, :);
        other = setdiff(1:num_rows, sampled);
        other_data = data(other,:);
        num_samples = numel(sampled);
        permed_index = [sampled, other];
    end
end
clear data;

% calculate the weights distances 
other_points = num_rows - num_samples;

if strcmp(opt.Metric, 'Euclidean')
    A = sqdist(sample_data',sample_data'); 
    B = sqdist(other_data',sample_data'); 
    
    if(opt.tau == -1) % automatic
        auto_tau = prctile(A(:),18); 
        tau = auto_tau*1.1;
        A = exp(-A/tau);
        B = exp(-B/tau);
    end
    if (opt.tau == -2) %z-p
        
        inv_d = zeros(num_rows,1);
        for i = 1:num_samples
            temp = [A(i,:) B(i,:)];
            [~,ind] = sort(temp,'ascend');
            inv_d(i) = 1/sqrt(temp(ind(opt.K)));
        end
        for i = 1:num_rows- num_samples
            temp = B(:,i);
            [~,ind] = sort(temp,'ascend');
            inv_d(i+num_samples) = 1/sqrt(temp(ind(opt.K)));
        end            
        tau_A = inv_d(1:num_samples)*inv_d(1:num_samples)';
        tau_B = inv_d(num_samples+1:end) * inv_d(1:num_samples)';
        A = exp(-A.*tau_A);
        B = exp(-B.*tau_B');
    end
    if opt.tau > 0
        A = exp(-A/tau);
        B = exp(-B/tau);
    end
    %A = A - diag(diag(A));
    %A(1:size(A,1)+1:size(A,1)*size(A,2)) = 1; %set diagonal to 1 
end

clear sample_data other_data;


% Normalize A and B using row sums of W, where W = [A B; B' B'*A^-1*B].
% Let d1 = [A B]*1, d2 = [B' B'*A^-1*B]*1, dhat = sqrt(1./[d1; d2]).

if opt.Laplacian == 'n'
%    B_T = B';
    d1 = sum(A, 2) + sum(B, 2);
    d2 = sum(B, 1)' + B'*(pinv(A)*sum(B, 2));
    dhat = sqrt(1./[d1; d2]);
    A = A .* (dhat(1:num_samples)*dhat(1:num_samples)');
    %B1 = dhat(1:num_samples)*dhat(num_samples+(1:other_points))';
    B = B .* (dhat(1:num_samples)*dhat(num_samples+(1:other_points))');
    clear d1 d2 dhat;
end


% Do orthogalization and eigendecomposition
%Asi = sqrtm(pinv(A));
[Bx, G] = eig(A);
G = diag(G);
Asi = Bx * diag(1./sqrt(G))*Bx';
%B_T = B';
BBT = B*B';
%W = double(zeros(size(A, 1)+size(B_T, 1), size(A, 2)));
%W(1:size(A, 1), :) = A;
%W(size(A, 1)+1:size(W, 1), :) = B_T;
W = [Bx * diag(sqrt(G)) ; B' * Bx * diag(1./sqrt(G))];
clear B;
% Calculate R = A + A^-1/2*B*B'*A^-1/2
R = A + Asi*BBT*Asi;
R = (R + R')/2; % Make sure R is symmetric, sometimes R can be non-symmetric because of numerical inaccuracy
[U,E] = eig(R);
[~, ind] = sort(diag(E), 'descend');
U = U(:, ind); % in decreasing order
E = E(ind, ind); % in decreasing order
E = diag(E);
%clear A R BBT;
%W = W*Asi;
%phi = bsxfun(@rdivide, W*U, sqrt(E'));
phi = W * Bx' * (U * diag(1./sqrt(E)));
phi(permed_index,:) = phi;
%phi = real(phi);
E = 1-E;

phi = phi(:,1:neig); 
E = E(1:neig);
end



function opt = check_opt(opt)
if(~isfield(opt, 'tau'))
    opt.tau = 40; 
end
if(~isfield(opt,'Laplacian'))
    opt.Laplacian = 'n'; 
end
if(~isfield(opt,'Metric'))
    opt.Metric= 'Euclidean'; 
end
if(~isfield(opt,'numsample'))
    opt.numsample = 500; 
end
if(~isfield(opt,'neig'))
    opt.neig = 500; 
end
end




