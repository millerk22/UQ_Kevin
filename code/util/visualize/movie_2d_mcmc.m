function movie_2d_mcmc(iter_stats)
%--------------------------------------------------------------------------
%  Description : 
%       Play a movie of the sampled function displayed on the Eigenvector
%       Coefficients of u
%
%  Input : 
%       u_hat_rec : record of u_hat, i.e., coefficient on eigenvectors
%       u_hat_mean_rec : record of mean of u_hat
%       Neig : number of components to display
%
%  Output : 
% -------------------------------------------------------------------------
u_hat_rec = iter_stats.u_hat_rec; 
u_hat_mean_rec = iter_stats.u_hat_mean_rec; 
max_iter = size(u_hat_rec, 3); 
clow = prctile(u_hat_rec(:), 0.01); 
chigh = max(u_hat_rec(:)) * 0.6; 


NClass = size(u_hat_rec, 2);
NEig = size(u_hat_rec, 1); 

H = figure;
set(H,'DoubleBuffer','on');
frame_rate = 40; 

for k=1:frame_rate:max_iter
    subplot(211); 
    imagesc(u_hat_rec(:, :, k), [clow, chigh]); 
    subplot(212); 
    imagesc(u_hat_mean_rec(:, :, k)); 
    title(['Iter : ', num2str(k)]); drawnow;
    clf; 
end

    subplot(211); 
    imagesc(u_hat_rec(:, :, k), [clow, chigh]); 
    subplot(212);
    imagesc(u_hat_mean_rec(:, :, k)); drawnow; 
    title(['Iter : ', num2str(k)]); 
    hold off;
end



