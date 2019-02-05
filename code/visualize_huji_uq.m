clear all
data_root = './data/HUJI/';         %Location of unshuffled data
temp_root = './data/HUJI/temp/';    % Shuffled data
load('iter_UQ_results_huji_181105_5.mat');
sol1  = class_result{2};
m1 = posterior_mean{2};
load([data_root,'ground_truth_huji_disney_60seg.mat'])
load([temp_root, 'train_test_index.mat'], 'train_index', 'test_index');
shuffle_index = [train_index test_index];
ground_truth = ground_truth(shuffle_index);
ground_truth    = ground_truth(1:5:length(ground_truth));
%%
%idx = 18000:36421;
idx = 1:size(ground_truth,1);
temp2 = sol1(idx);
gt = ground_truth(idx);
uncertainty = var(m1(idx,:)');

[~,idx] = sort(uncertainty','descend');
gt = gt(idx);
temp2 = temp2(idx);
%%
acc = cumsum(temp2 == gt)./(1:size(temp2,1))';
plot((1:size(temp2,1)) * 100/size(temp2,1),acc*100, 'linewidth', 2);
xticks([0:10:100])
yticks([0:5:100])
ylim([76,100])
xlim([0,100])
xtickformat('percentage');
ytickformat('percentage');

set(gca, 'fontsize', 12,'Fontname','Times')
xlabel('Confidence percentile','fontsize',14,'Interpreter','latex');
ylabel('Classification accuracy','fontsize',14,'Interpreter','latex');
%% Separate interval
clear acc perc
for i = 1:20
    inc = floor(size(ground_truth,1)/20);
    acc(i) = sum(temp2(1 + (i-1) * inc: i*inc) == gt(1 + (i-1) * inc: i*inc))/inc;
    perc(i) = (i-0.5) * inc / size(ground_truth,1);
end
%perc = [perc 1];
%acc = [acc acc(end)];
%perc = [0 perc];
%acc = [1  acc];
%acc(i+1) = sum(temp2(1 + i *  inc:end) == gt (1 + i *  inc:end)) / (36421 - i * inc);
%perc(i+1) = 1;
%stairs(perc * 100,acc,'-', 'linewidth', 2)
bar(perc*100,acc*100,1);
axis tight

xticks([0:20:100])
yticks([0:20:100])
ylim([0,100])
xlim([0,100])
xtickformat('percentage');
ytickformat('percentage');
set(gca, 'fontsize', 18,'Fontname','Times')
xlabel('Top $x\%$ confident data points','fontsize',24,'Interpreter','latex');
ylabel('Classification accuracy','fontsize',24,'Interpreter','latex');
%%
correct = find(temp2 == gt);
uncorrect = find(temp2~= gt);
figure 
[a,b] = hist(uncertainty(correct),15);
a = a / numel(correct);
bar(b,100 * a,1)
xticks([0:0.03:0.15])
%ytickformat('percentage');
axis tight
set(gca, 'fontsize', 12,'Fontname','Times')
xlabel('Confidence score, $S(i)$','fontsize',14,'Interpreter','latex');
ylabel('Percentage, \%','fontsize',14,'Interpreter','latex');
title('Correct classification', 'fontsize',14,'Interpreter','latex')
%%
figure 
[a,b] = hist(uncertainty(uncorrect),15);
a = a / numel(uncorrect);
bar(b,100 * a,1)
xticks([0:0.03:0.15])
%ytickformat('percentage');
axis tight
set(gca, 'fontsize', 12,'Fontname','Times')
xlabel('Confidence score, $S(i)$','fontsize',14,'Interpreter','latex');
ylabel('Percentage, \%','fontsize',14,'Interpreter','latex');
title('Wrong classification', 'fontsize',14,'Interpreter','latex')
