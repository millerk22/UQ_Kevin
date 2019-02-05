function fid_list = get_fid_list(fid)
% Return a list of fidelity points
% fid{i} = list of fidelity points of class i
    fid_list = [];
    for i =1:numel(fid)
        fid_list = [fid_list; fid{i}];
    end
end