%% Graph Construction
rng(2)
digits      = [1, 4, 7, 9];
NClass      = numel(digits); 
num_points  = [500, 500, 500,500];
disp(['Generating Graph']);
tic;
Neig        = 300;
ylbs        = 1:NClass;
[V, E, outargs]     = makemnistgraph(digits, num_points, Neig, 's');
time_       = toc;
disp(['Generating graph took ', num2str(time_),  ' s']);
im          = outargs.images;
ground_truth        = outargs.ground_truth;
%%
% figure;
% plot(E);
% title('Eigenvalues of graph Laplacian');
% drawnow;
% % plot to visualize
% cl_dots = 'rbgkc';
% figure;
% plot(E);
% figure;
% for i = 1:numel(ylbs)
%     scatter(V(ground_truth == ylbs(i), 2), V(ground_truth == ylbs(i), 3), [cl_dots(i), '.']) ;
%     hold on;
% end
% hold off;
% title('Projection onto 2nd and 3rd eigenvector');
% drawnow;
% 
% [y_sp, ~] = kmeans(V(:, 2:3), 3, 'Distance', 'cosine'); 
% for i = 1:numel(ylbs)
%     scatter(V(y_sp == ylbs(i), 2), V(y_sp == ylbs(i), 3), [cl_dots(i), '.']) ;
%     hold on;
% end
% acc_sp = unsupervised_accuracy(y_sp, ground_truth); 
% title(['Spectral Clustering Accuracy = ', num2str(acc_sp)]);
% drawnow; 


%% Running MCMC
fidelity_percent = 0.03; % generate fidelity
N = sum(num_points); 
Ntemp = 0; 
fid = {}; 
for i = 1:NClass
    class_i = find(ground_truth == i);
    class_i = class_i(randperm(length(class_i)));
    fid{i} = class_i(1:ceil(fidelity_percent * length(class_i)));
end
disp(['Number of points: ', num2str(N)]);
disp(['Percent of fidelity: ', num2str(fidelity_percent)]);


disp('Starting MCMC...');
tic;
gamma=0.1; % obs noise std
max_iter=2*10^4; % number of mcmc steps
%Phi_eval = get_multiclass_likelihood_function('bls', gamma);
%E(1) = 0.2; 
[m1, iter_stats1] = mcmc_gibbs_probit_multiclass(max_iter, gamma,  ...
                            V, E,fid); 
[m2, iter_stats2] = mbo_multiclass(0.1, 1, 50, ...
                            V, E, fid); 
%[m1,iter_stats1] = mcmc_pcn_core_multiclass(beta, max_iter, Phi_eval, V,E,fid);
time_ = toc;
disp(['Running MCMC took ', num2str(time_),  ' s']);
%disp(['Average Acceptance Probability = ', num2str(iter_stats1.mcmc_cum_prob(end))]);

%%
[~, yp] = max(m1, [], 2); 
[~, ypMBO] = max(m2, [], 2);
to_eval = 1:N;
for i = 1:NClass
    to_eval = setdiff(to_eval, fid{i});
end
disp(['Classification Accuracy = ', num2str(sum(yp(to_eval) == ground_truth(to_eval))/N)]);

disp(['Classification Accuracy MBO= ', num2str(sum(ypMBO(to_eval) == ground_truth(to_eval))/N)]);
%%

confidence = var(m1, [],2);
highCandidate    = {[135,   151,    183,    232],
                    [603,   613,    686,    774],
                    [1125,  1166,   1200,   1394],
                    [1622,  1699,   1729,   1975]
                    };
lowCandidate    =  {[224,   439,    345,    305],
                    [680,   914,    835,    917],
                    [1163,  1205,   1319,   1321],
                    [1598,  1758,   1797,   1834]
                    };
fig=  figure;
subplot = @(m,n,p) subtightplot (m, n, p, [0.07 -0.0], [0.01 0.05], [0 0]);
for class = 1:NClass
    for i = 1:4
        subplot(4,4, (class -1) * 4 + i)
        img = reshape(im(highCandidate{class}(i),:),[28,28]);
        imshow(img,[], 'InitialMagnification', 500);
        title(sprintf('$S(i)=$ %.3f', confidence(highCandidate{class}(i))),...
             'fontsize',12, 'interpreter', 'latex')
    end
end
%%

%draw pics with lower 5% confidence
confidence = var(m1, [],2);
for class = 1:NClass
    idx_per_class = find(ground_truth == class);
    [~,idx] = sort(confidence(idx_per_class), 'descend');
    for i  = 1:10
        img = reshape(im(idx_per_class(idx(i)),:),[28,28]);
        figure('rend','painters','pos',[10 10 30 30])
        imshow(img,[], 'InitialMagnification', 400);
        title(sprintf('$S(i)=$ %.3f', confidence(idx_per_class(idx(i)))),...
             'fontsize',14, 'interpreter', 'latex')
        %title(['$S(i)=$', num2str(confidence(idx_per_class(idx(i))))],'interpreter', 'latex');
        saveas(gcf,['./visualization/MNIST/higher_confidence/digit_', num2str(digits(class)), '_image_', ...
            int2str(idx_per_class(idx(i))),'.fig']);
        close all
    end
    [~,idx] = sort(confidence(idx_per_class), 'ascend');
    for i  = 1:10
        img = reshape(im(idx_per_class(idx(i)),:),[28,28]);
        figure('rend','painters','pos',[10 10 30 30]);
        imshow(img,[], 'InitialMagnification', 400);
         title(sprintf('$S(i)=$ %.3f', confidence(idx_per_class(idx(i)))),...
             'fontsize',14, 'interpreter', 'latex')
        saveas(gcf,['./visualization/MNIST/lower_confidence/digit_', num2str(digits(class)), '_image_', ...
            int2str(idx_per_class(idx(i))),'.fig']);
        close all
    end
end
%%
[~, idx] = sort(confidence, 'ascend');
percent = 0.05
num = N * percent
for i=1 : num
    img = reshape(im(idx(i),:),[28,28]);
    
    figure;
    imshow(img,[]);
    title(['$S(i)=$', num2str(confidence(idx(i)))], 'interpreter', 'latex');
    saveas(gcf,['./visualization/MNIST/lower_confidence/image_',int2str(idx(i)),'.jpg']);
    close all
end 
[~, idx] = sort(confidence, 'descend');
for i=1 : num
    img = reshape(im(idx(i),:),[28,28]);
    figure;
    imshow(img,[]);
    title(['$S(i)=$', num2str(confidence(idx(i)))], 'interpreter', 'latex');
    saveas(gcf,['./visualization/MNIST/higher_confidence/image_',int2str(idx(i)),'.jpg']);
    close all
end 


%%

confidence = var(m1, [],2);
correct = find(yp == ground_truth);
for i = 1:NClass
    correct = setdiff(correct, fid{i});
end
uncorrect = find(yp~=ground_truth);
figure 
[a,b] = hist(confidence(correct),15);
a = a / numel(correct);
bar(b,100 * a,1)
xticks([0:0.03:0.2])
%ytickformat('percentage');
axis tight
set(gca, 'fontsize', 20,'Fontname','Times')
xlabel('Confidence score, $S(i)$','fontsize',24,'Interpreter','latex');
ylabel('Percentage, \%','fontsize',24,'Interpreter','latex');
title('Correct classification', 'fontsize',24,'Interpreter','latex')
%%
figure 
[a,b] = hist(confidence(uncorrect),15);

a = a / numel(uncorrect);
bar(b,100 * a,1)
xticks([0:0.03:0.2])
%ytickformat('percentage');
axis tight
set(gca, 'fontsize', 20,'Fontname','Times')
xlabel('Confidence score, $S(i)$','fontsize',24,'Interpreter','latex');
ylabel('Percentage, \%','fontsize',24,'Interpreter','latex');
title('Wrong classification', 'fontsize',24,'Interpreter','latex')


