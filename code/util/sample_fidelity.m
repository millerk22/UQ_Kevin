function [fid] = sample_fidelity(gt, training_size, fid_perc, varargin)
type = unique(gt); % 1-13
N = length(gt);
opts = check_opt(varargin{1});
for i = 1:length(type)
    idx = find(gt == type(i));
    idx = idx(idx <= training_size * N);
    m = length(idx);
    switch opts.type
        case 'random'
            idx = idx(randperm(m));
            idx = idx(1:ceil(fid_perc * m));
        case 'uniform'
            stride = floor ( 1 / fid_perc);
            idx = idx(1:stride:m);
    end
    fid{i} = idx;
end
end

function opt = check_opt(opt)
    if(~isfield(opt, 'type')); 
        opt.type = 'random'; 
    end
end