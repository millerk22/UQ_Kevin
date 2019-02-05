function movie_nodes_mcmc(iter_stats)
%--------------------------------------------------------------------------
%  Description : 
%       Play a movie of the sampled function displayed on the Nodes of 
%       the graph. 
%
%  Input iter_stats with field: 
%       y_mean_rec : record of the mean of the y-prediction. 
%       u_rec : record of the current sample of u. 
%
%  Output : 
% -------------------------------------------------------------------------
y_mean_rec = iter_stats.y_mean_rec; 
u_rec = iter_stats.u_rec; 
N = size(y_mean_rec,1); 
max_iter = size(y_mean_rec, 2); 
w1=[0,N,-1,1]; 
w2=[0,N,-1.0,1.0]; 
xx=1:N;
H = figure;
set(H,'DoubleBuffer','on');
frame_rate = 20; 
for k=1:frame_rate:max_iter
    subplot(211), axis(w1); hold on;
    plot(xx,u_rec(:,k),'ro');
    subplot(212), axis(w2); hold on;
    plot(xx,y_mean_rec(:,k)','+'); drawnow; 
    clf; hold off;
end

subplot(211), axis(w1); hold on;
plot(xx,u_rec(:,k)','ro');
subplot(212), axis(w2); hold on;
plot(xx,y_mean_rec(:,k)','+'); drawnow; hold off;

end


