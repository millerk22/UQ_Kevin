function y = Nind2vec(vec_ind, varargin)
% Convert an array of integer indices to its vector form
% e.g. [1, 2, 1] ->  [0, 1; 1, 0; 0, 1]; 
ind = unique(vec_ind);
if numel(varargin) == 1
    NClass = varargin{1}; 
else
    NClass = max(ind); 
end
N = numel(vec_ind);
vec_ind = reshape(vec_ind, 1, N); 
y = zeros(NClass, N);
indt = 1:NClass:N * NClass;
indt = indt + vec_ind - 1;
y(indt) = 1;
y = y';
end