% construct graphs for for MNIST using method in Hu et al paper
% pca the digits to 50 dimensions
% pick 10 nearest neighbors as graphs
% compute using Perona Scaling.

function [V, E, outargs] = makemnistgraph(digits, num_points, Neig, Ltype)
global initial_path
num_components = 50;
num_nbs = 15;
[lb, im] = mnist_data(digits);

% subsample a smaller portion, 40 % totally 4000 points
Ndigits = numel(digits); 
N = zeros(Ndigits, 1); 
n = num_points; 
for i = 1:Ndigits
    N(i) = sum(lb == digits(i));
    if(n(i) >= N(i))
        n(i) = N(i);
    end
end
ind = zeros(sum(n), 1); 
ntemp = 0; 
Ntemp = 0; 
for i = 1:Ndigits
    ind((1+ntemp):(n(i)+ntemp)) = randperm(N(i), n(i)) + Ntemp; 
    ntemp = ntemp + n(i); 
    Ntemp = Ntemp + N(i); 
end
lb = lb(ind);

% do PCA onto 50 dimensions first.
% sc is the projection of the digits onto the axis. (N x k) matrix
[~, sc] = pca(im, 'NumComponents', num_components);
sc = sc(ind, :);
im = im(ind, :); 

% do 10 nearest neighbors, and mean dist between 10 nearest neighbors
num_points = size(sc, 1);
dist = sqdist(sc', sc');
[dist, index] = sort(dist, 2, 'ascend');
d_sp = dist(:, 2:num_nbs+1);
j_sp = index(:, 2:num_nbs+1);
clear dist index;

% compute the weights via the scaling by mean of closest 10 dist
dsum_sp = sum(d_sp, 2);
dmean_sp = dsum_sp / num_nbs;
w_sp = bsxfun(@rdivide, d_sp, dmean_sp);
w_sp = exp(-(w_sp .* w_sp)/ 3);


if Ltype == 'u'
    % compute and store sparse matrix L(unnormalized)
    i_sp = reshape((1:num_points)' * ones(1, num_nbs), 1, num_points * num_nbs , 1);
    j_sp = reshape(j_sp, numel(j_sp), 1);
    W = sparse(i_sp, j_sp, w_sp);
    W = .5 * (W + W'); % symmetrize
    wsum_sp = sum(W, 2);
    wsum_sp(wsum_sp < 1e-6) = 1e-6;
    isum_sp = 1:num_points;
    D = sparse(isum_sp, isum_sp, wsum_sp);
    Dsqrtinv =  sparse(isum_sp, isum_sp, 1./sqrt(wsum_sp));
    L = D - W;
end

if Ltype == 's'
    % compute and store sparse matrix L(symmetric)
    i_sp = reshape((1:num_points)' * ones(1, num_nbs), 1, num_points * num_nbs , 1);
    j_sp = reshape(j_sp, numel(j_sp), 1);
    W = sparse(i_sp, j_sp, w_sp);
    W = .5 * (W + W'); % symmetrize
    wsum_sp = sum(W, 2);
    wsum_sp(wsum_sp < 1e-6) = 1e-6;
    isum_sp = 1:num_points;
    D = sparse(isum_sp, isum_sp, wsum_sp);
    Dsqrtinv =  sparse(isum_sp, isum_sp, 1./sqrt(wsum_sp));
    L = speye(num_points) - Dsqrtinv * W * Dsqrtinv;
end

% compute the eigenvectors
[V, E] = eigs(L, Neig, 'SM');
E = diag(E);

ground_truth = ones(size(lb));
if Ndigits == 2
    ground_truth(lb == digits(2)) = -1;
else
    for i = 1:Ndigits
        ground_truth(lb == digits(i)) = i; 
    end
end


% % plot to visualize
% cl_dots = 'rbgkc'; 
% figure;
% plot(E);
% figure;
% if Ndigits == 2
% plot(V(ground_truth > 0, 2), V(ground_truth > 0, 3), 'r.');
% hold on;
% plot(V(ground_truth < 0, 2), V(ground_truth < 0, 3), 'b.');
% hold off;
% else
%     for i = 1:Ndigits
%         scatter(V(ground_truth == i, 2), V(ground_truth == i, 3), [cl_dots(i), '.']) ;
%         hold on; 
%     end
%     hold off; 
% end
% 
% figure;
% plot(V(:, 3));

if Ltype == 's'
    tag = 'sym';
else
    tag = 'unnorm'; 
end

dstr = ''; 
for i = 1:Ndigits
    dstr = [dstr, num2str(digits(i))]; 
end
outargs = {}; 
outargs.images = im; 
outargs.ground_truth = ground_truth; 
% fname = [initial_path, '/data/MatFiles/mnist', dstr 'graph_', tag, '.mat']; 
% save(fname, 'im', 'j_sp', 'i_sp', 'w_sp', 'L', 'V', 'E', 'ground_truth');

end




