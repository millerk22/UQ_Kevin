function [ Phi_eval ] = get_multiclass_likelihood_function(type, varargin)
%--------------------------------------------------------------------------
% Author: Xiyang Luo <xylmath@gmail.com> , UCLA
%
% This file is part of the diffuse-interface graph algorithm code.
% There are currently no licenses.
%
%--------------------------------------------------------------------------
%  Description :
%       Binary Negative Log Likelihood function handle generator. Add all
%       new log likelihood functions to this file.
%  Input :
%`      type : str, name of the likelihood function. E.g., 'probit'
%       varargin : parameters for log likelihood.
%
%  Output :
%       Phi_eval : function handle, a function handle that maps the triplet
%           (u, fid, y) to a real number, i.e., the negative log
%           likelihood. u = N x 1 vector. fid is cell array with fid{1}
%           being indices of fidelity points for 1st class. y being the
%           label set. Default [-1, 1].
% -------------------------------------------------------------------------

if strcmp(type, 'bls')
    gamma = varargin{1};
    Phi_eval = @(u, fid, y) bls_loglike(u, fid, y, gamma);
end
if strcmp(type, 'bls_normalize')
    gamma = varargin{1};
    Phi_eval = @(u, fid, y) bls_loglike_normalize(u, fid, y, gamma);
end
if strcmp(type, 'logit')
    gamma = varargin{1};
    Phi_eval = @(u, fid, y) logit_loglike(u, fid, y, gamma);
end
if strcmp(type, 'Ginzberg')
    gamma = varargin{1};
    eta = varargin{2};
    Phi_eval = @(u, fid, y) Ginzberg(u, fid, y, gamma, eta);
end
end


function res = bls_loglike(u, fid, label, gamma)
res = 0; 
for j = 1:numel(fid)
    [~, temp] = max(u(fid{j}, :), [], 2);
    res = res + sum(temp ~= label(j))*sqrt(2)/(2*gamma^2);
end
end

function res = bls_loglike_normalize(u, fid, label, gamma)
res = 0; 
for j = 1:numel(fid)
    [~, temp] = max(u(fid{j}, :), [], 2);
    res = res + mean(temp ~= label(j))*sqrt(2)/(2*gamma^2);
end
end

function res = logit_loglike(u, fid, y, gamma)
res = 0;
for i = 1:numel(fid)
    ind = fid{i};
    for j = 1:numel(ind)
        p1 = sum(exp(-u(ind(j), :) / gamma)); 
        p2 = exp(-u(ind(j), i) / gamma); 
        res = res + log(p2/p1);
    end
end
end

function res = Ginzberg(u, fid, y, gamma, eta)
res = 0;
u = u + max(-min(min(u)),0);
rsum = sum(u,2);
u = bsxfun(@rdivide, u, rsum);
for i = 1:numel(fid)% fidelity term
    ind = fid{i};
    for j = 1:numel(ind)
        temp = u(ind(j),:);
       temp(:, i) = temp(:, i)-1;
        res = res + sum(sum(temp.^2)/gamma);
    end
end

    mw = ones(size(u,1),1);
    for j = 1:numel(y)% multiwell potential
        u(:,j) = u(:,j) - 1;
        mw = mw.* (((sum(abs(u),2)).^2)/4);
        %res = res + sum(sum(abs(u))/eta);
        u(:, j)=u(:,j)+1;
    end 
    res =res + sum(mw)*eta;

end


