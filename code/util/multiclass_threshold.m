function v = multiclass_threshold(u)
%%
% u(u<0) = 0;
% u(u>1) = 1;
% for i =1:size(u,2)
%     u(:,i) = u(:,i)/sum(u(:,i));
% end
%%
v = zeros(size(u'));
[~, ind] = max(u, [], 2);
ind = ind + size(u,2)*(0:1:(size(u,1)-1))';
size(u);
v(ind) = 1;
v = v';
end
