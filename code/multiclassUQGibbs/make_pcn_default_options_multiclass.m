function opt = make_pcn_default_options_multiclass(opt)
%--------------------------------------------------------------------------
% Author: Xiyang Luo <xylmath@gmail.com> , UCLA 
%
% This file is part of the diffuse-interface graph algorithm code. 
% There are currently no licenses. 
%
%--------------------------------------------------------------------------
%  Description :
%       Fill in missing values in options with default value for pCN
%       algorithm. 
%
%       Options Fields: 
%           isrec_u: bool, when true, u is recorded in output. 
%                           Default True.
%           Verbose: bool, when true, displays iteration number along with
%                           acceptance probability. Default True. 
%           ZeroMean: bool, when true, u is constrained to satisfy mean(u)
%                           = 0. Default True. 
%           RandomBatchSize: int, Batch size of the standard Gaussian
%                           variable z_k used. Larger Batch is faster but
%                           more memory intensive. Default -1, i.e., size
%                           of batch is the entire sequence length MaxIter.
%           LambdaBar: float/str. If Inf, spectral projection is used.
%                       i.e.,  u is projected to space spanned by V. If <
%                       Inf, approximate the rest of the eigenvalues by
%                       LambdaBar. If 'auto', then LambdaBar is set to be
%                       max(E), i.e, largest of current eigenvalue. 
%           BurnIn: int, Number of samples to discard during initial phase 
%                   of MCMC
%           VerboseWindowLength: int, Length of interval to report summary
%                               statistics for 'Verbose' option
%           Init: str, modes of initialization. Either 'fid', or 'rand'.
%                   'fid' initializes by fidelity. 
% -------------------------------------------------------------------------
if(~isfield(opt,'isrec_u'))
    opt.isrec_u = false;
end
if(~isfield(opt,'Verbose'))
    opt.Verbose = false;
end
if(~isfield(opt,'ZeroMean'))
    opt.ZeroMean = false;
end
if(~isfield(opt,'RandomBatchSize'))
    opt.RandomBatchSize = -1; 
end
if(~isfield(opt,'T'))
    opt.T = 1.0; 
end
if(~isfield(opt,'LambdaBar'))
    opt.LambdaBar = Inf; 
end
if(~isfield(opt,'BurnIn'))
    opt.BurnIn = 0; 
end
if(~isfield(opt,'VerboseWindowLength'))
    opt.VerboseWindowLength= 500; 
end
if(~isfield(opt,'Init'))
    opt.Init= 'fid'; 
end
if(~isfield(opt,'AdditionalAvgFunc'))
    opt.AdditionalAvgFunc = {}; 
end

end