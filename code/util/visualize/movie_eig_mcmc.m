function movie_eig_mcmc(u_hat_rec,u_hat_mean_rec,Neig)
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
max_iter = size(u_hat_rec, 2); 
w1=[0,Neig,-1,1]; 
w2=[0,Neig,-1,1]; 
xx=1:Neig;
H = figure;
set(H,'DoubleBuffer','on');
frame_rate = 20; 

for k=1:frame_rate:max_iter
    subplot(211), axis(w1); hold on;
    plot(xx,sign(u_hat_rec(1:Neig,k))','ro');
    subplot(212), axis(w2); hold on;
    plot(xx,u_hat_mean_rec(1:Neig,k)','+'); drawnow; 
    clf; hold off;
end

    subplot(211), axis(w1); hold on;
    plot(xx,sign(u_hat_rec(1:Neig,k))','ro');
    subplot(212), axis(w2); hold on;
    plot(xx,u_hat_mean_rec(1:Neig,k)','+'); drawnow; 
    hold off;
end


