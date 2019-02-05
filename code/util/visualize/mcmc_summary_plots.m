function mcmc_summary_plots(ratio, ct, u_hat_rec, u_hat_mean_rec, ...
    u_rec,y_mean_rec, pointwise_err_rec, Neig)
%--------------------------------------------------------------------------
%  Description : 
%       Make a summary plot of all iteration statistics from the MCMC
%       iterations.
%
%  Input : 
%       y_mean_rec : record of the mean of the y-prediction. 
%       u_rec : record of the current sample of u. 
%       pointwise_err_rec : record of mean prediction error from sample.
%       Neig : number of eigenvector components to display
%       u_hat_rec : record of u_hat, i.e., coefficient on eigenvectors
%       u_hat_mean_rec : record of mean of u_hat
%       ratio : overall classification rate per timestep
%       ct : acceptance probability of the MCMC. 
%
%  Output : 
% -------------------------------------------------------------------------

N = size(y_mean_rec,1); 
max_iter = size(y_mean_rec,2); 

figure;
title('accuracy vs time');
plot(ratio);
figure;
title('acceptance rate');
plot(ct);
xx=1:N;
w1=[0,20,-12.5,12.5]; w2=[0,20,-1,1]; w3=[0,20,0.0,12.5];

figure;
subplot(211), axis(w1); hold on;
plot(xx(1:20),u_hat_rec(1:Neig,max_iter)','ro');
subplot(212), axis(w2); hold on;
plot(xx(1:20),u_hat_mean_rec(1:Neig,max_iter)','+');
hold off;
%    subplot(313), axis(w3); hold on;
%    plot(xx(1:20),et(1:20,max_iter)','g*'); drawnow; hold off;

w1=[0,N,-1.,1.]; w2=[0,N,-1.0,1.0]; w3=[0,N,0.0,1.0];

figure;
subplot(311), axis(w1); hold on;
plot(xx,sign(u_rec(:,max_iter))','ro');
subplot(312), axis(w2); hold on;
plot(xx,y_mean_rec(:,max_iter)','+');
subplot(313), axis(w3); hold on;
plot(xx,pointwise_err_rec(:,max_iter)','g*'); drawnow; hold off;
end