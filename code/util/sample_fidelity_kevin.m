function [fid] = sample_fidelity_kevin(gt, training_size, fid_perc, varargin)
type = unique(gt); % the different classes
N = length(gt);
opts = check_opt(varargin{1}); % check if have the opt.type variable

for i = 1:length(type)          %for each class
    idx = find(gt == type(i));  % find indices of the current class
    idx = idx(idx <= training_size * N); % restrict to indices that are in "training chunk" 
    m = length(idx);            % # of training indices in this class
    switch opts.type
        case 'random'           % if random, take random subset of these 
            idx = idx(randperm(m));
            idx = idx(1:ceil(fid_perc * m));
        case 'uniform'          % if uniform take every 'stride'-th element
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