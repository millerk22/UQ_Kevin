function fig = tile(X,varargin)
if numel(varargin) == 1
    opt = varargin{1}; 
else
    opt = {}; 
end
opt = check_opt(opt); 
fig = figure; 
p = size(X,ndims(X)); 
m = ceil(sqrt(p)); 
n = ceil(p/m); 
for i = 1:p
    subplot(m,n,i); 
    if(ndims(X) == 2)
        v = reshape(X(:,i),opt.xsize, opt.ysize); 
    else
        v = X(:,:,i); 
    end
    imagesc(v); 
    axis off; 
    switch opt.colormap
        case 'gray'
            colormap gray; 
    end
end
end

function opt = check_opt(opt)
if(~isfield(opt,'colormap'))
    opt.colormap = 'default';
end
if(~isfield(opt,'num_rows'))
    opt.num_rows = 'auto'; 
end
end

