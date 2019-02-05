function outstr =  struct_description(struct)
    outstr = {};
    fn = fieldnames(struct);
    for i = 1:numel(fn)
        outstr{i} = [fn{i},'_',num2str(getfield(struct, fn{i})),'_'];
    end
    outstr = strjoin(outstr);
end