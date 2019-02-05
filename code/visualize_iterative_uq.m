% Visualize iterative UQ
function visualize_iterative_uq(suffix, field)
iter_filename       = ['iter_UQ_results_huji', suffix, '.mat'];
noniter_filename    = ['noniter_UQ_results_huji', suffix, '.mat'];
% class_result(_noniter):     niter cell, classification from UQ
% class_result(_noniter)_MBO: niter cell, classification from MBO
% fid_perc:                   niter+1 cell,
% posterior_mean(_noniter)    niter cell, 
% summary_stats(_noniter):    niter cell, stats of classification from UQ,
% summary_stats(_noniter)_MBO:niter cell, stats of classification from MBO
% summary_stats: 
%   -C:                         unnormalized confusion matrix
%                               ground_truth x classification
%   -recall/precision_confmat:  confusion matrix normalized accordingly
%   -recall:                    mean recall over seven classes on the test
%                               set
%   -precision                  mean precision over seven classes on the
%                               test set
%   -acc                        accuracy on the test set
load(noniter_filename)
load(iter_filename)

niter       = numel(class_result);
fid_perc    = fid_perc(1:niter); % remove the last trailing fid_perc
%% Accuracy vs fidelity percentage
fig         = figure();
hold on
if strcmp(field, 'recall')
    get_stats   =  @(x) x.recall;
elseif strcmp(field, 'acc')
    get_stats   =  @(x) x.acc;
else
    get_stats   =  @(x) x.precision;
end
linestyle   =  {'--x','--*','-o','-+'};
linewidth   =   1.5;
plt(1)      =   plot(fid_perc * 100, cellfun(get_stats, summary_stats_MBO) * 100, ...
                linestyle{1}, 'linewidth', linewidth);
plt(2)      =   plot(fid_perc* 100, cellfun(get_stats, summary_stats)* 100, ...
                 linestyle{2}, 'linewidth', linewidth);
plt(3)      =   plot(fid_perc* 100, cellfun(get_stats, summary_stats_noniter_MBO)* 100, ...
                linestyle{3}, 'linewidth', linewidth);
plt(4)      =   plot(fid_perc* 100, cellfun(get_stats, summary_stats_noniter)* 100,...
               linestyle{4}, 'linewidth', linewidth);
axis tight
xticks([6:6:30])
legend({'MBO-iter','UQ-iter','MBO', 'UQ'})
set(gca, 'fontsize', 18,'Fontname','Times')
xlabel('Fidelity Percentage, \%','interpreter','latex', 'fontsize',24)
ylabel('Accuracy, \%','interpreter','latex', 'fontsize',24)
end
