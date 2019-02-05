function [ u, u_hat] = mbo_multiclass(dt, eta, max_iter, V, E, fid, varargin)
s = 10; % inner loop size
N = size(V, 1);

NClass = numel(fid);
dt2 = dt / s;
mul = ones(size(E)) + dt2 * E;
mul = 1./mul;
%mul = ones(size(E)) - dt2 * E;
%MBO initialization
%randind = randi(NClass, N, 1);
%u = Nind2vec(randind, NClass);
u = zeros(N, NClass);

% -------- Original Code --------
for i = 1:NClass
    u(fid{i}, :) = 0;
    u(fid{i}, i) = 1;
end

%%

du_J = zeros(N, NClass);

for i = 1:max_iter
    for j = 1:s
        for k = 1:NClass
            du_J(fid{k}, :) = u(fid{k}, :);
            du_J(fid{k}, k) = du_J(fid{k}, k) - 1; 
        end
        d_hat = V' * du_J;
        u_hat = V' * u;
        [~,w] = size(u_hat);
        vec = zeros(size(u_hat));
        for index = 1 : w
            vec(:,index) = u_hat(:,index) .* mul;
        end
        
        u_hat = vec - dt2 * eta * d_hat;
        u = V*u_hat; 
    end
    u = multiclass_threshold(u); 
    for i = 1:NClass
        u(fid{i}, :) = 0;
        u(fid{i}, i) = 1;
    end
end

end

