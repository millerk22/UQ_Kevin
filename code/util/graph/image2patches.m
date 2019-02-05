function [ res ] = image2patches( image, patch_size,opt )

%--------------------------------------------------------------------------
% Author: Xiyang Luo <xylmath@gmail.com> , UCLA 
%
% This file is part of the diffuse-interface graph algorithm code. 
% There are currently no licenses. 
%
%--------------------------------------------------------------------------
% Description: Returns all patches extracted form image
%
% Inputs: 
%       image : Input Image, either gray or color
%       patch_size : Integer, half the size of the actual patch. nsz =
%                   2*patch_size+1
%       opt : Options
%           output_dim : output dimension, 2,3,or 4. 
%           kernel : flag for whether to apply spatial mask to patch. true
%                    or false. 
%           kernel_params : any parameters for making a kernel. 
%
% Outputs:
%       res : a 2D or 3D or 4D vector containing all image patches. 
%               (spatial, patch_dimension, channel)
%               or (spatial, patch*channel)
%
%--------------------------------------------------------------------------


if(max(image(:))>1)
    image = image/max(image(:)); 
end
opt = check_opt(opt); 
% convert an image to an array of patches
if(ndims(image) == 2); 
    image = reshape(image,[size(image,1), size(image,2),1]); 
end
im = padarray(image,[patch_size,patch_size],'symmetric','both'); % symmetric reflection
nsz = 2*patch_size+1; 
res = zeros(size(image,1),size(image,2),nsz*nsz,size(image,3)); 
for i = 1:size(image,1)
    for j = 1:size(image,2)
        res(i,j,:,:) = flatten_2d(im(i:i+2*patch_size, j: j+2*patch_size,:)); 
    end
end

if(opt.kernel) % multiply a mask on top of each patch
    mask = make_kernel(opt.kernel_params,patch_size); %2d array
    res = bsxfun(@times, res, reshape(mask(:), 1,1,[],1)); 
end

if(opt.output_dim == 3); 
    %flatten the result array res into 3-D
    res = reshape(res,size(res,1)*size(res,2),size(res,3),[]);    
end
if(opt.output_dim == 2); 
    %flatten the result array res into 2-D
    res = reshape(res,size(res,1)*size(res,2),size(res,3),[]);  
    res = permute(res,[2,3,1]); 
    res = reshape(res,size(res,1)*size(res,2),size(res,3)); 
    res = permute(res,[2,1]); 
end

end


function res = flatten_2d(arr)
    res = reshape(arr,size(arr,1)*size(arr,2),size(arr,3)); 
end


function opt = check_opt(opt)
    if(~isfield(opt,'output_dim'))
        opt.output_dim = 2; 
    end
    if(~isfield(opt,'kernel'))
        opt.kernel = false; 
    else
        if(~isfield(opt,'kernel_params'))
            opt.kernel_params = {}; 
        end
    end
end



function res = make_kernel(kernel_params,patch_size)
% make a 2d-kernel
nsz = patch_size*2+1; 
res = zeros(nsz, nsz); 
for i = 1:nsz
    for j = 1:nsz
        d = abs(i-patch_size-1)+abs(j-patch_size-1); 
        res(i,j) = 1/(1+2*d/patch_size); 
    end
end
res = sqrt(res); 
res = res./max(abs(res(:))); 
end

