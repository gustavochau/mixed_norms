clc;
clear;
close all;

N = 20;
num_groups = 50;
y = rand(N,1);

% create indices
% rng(6)
for ii=1:num_groups
    p = randperm(N);
    g_index{ii} = p(1:randi(N/2+1,1,1));
end

rho=1;

G = zeros(num_groups,N);
for ii=1:num_groups
    G(ii,g_index{ii})=1;
end

G_aux = [];
for ii=1:num_groups
%     G_aux = [G_aux; spdiag(G(ii,:))];
    G_aux = [G_aux; spdiags(G(ii,:)', 0, N, N)];
end    
% G_aux = sparse(G_aux);

lambda = 0.1;
params.y = y;
params.mat_form = @(w) reshape(w,[num_groups,N]);
    params.G_aux = G_aux;
    params.lambda = lambda;
    

%%%%%%%%%%%
%% x update
%%%%%%%%%%%

% z = rand(num_groups,N);
% v = abs(rand(num_groups,N));
% 
% 
% cvx_begin
%     variable x_cvx(N,1)
%     minimize( 0.5*sum_square( y-x_cvx) + 0.5*rho*sum_square(vec(z-G*diag(x_cvx)+v)))
% cvx_end
% 
% A = @(W) G*W;
% AT = @(W) G'*W;
% 
% z=vec(z);
% v=vec(v);
% 
% [ x_admm ] = xupdate_1inf_overlap( y,z,A,1,1,v,rho,params, AT );
% 
% max(abs(x_admm-x_cvx))

% %%%%%%%%%%%
% %% z update
% %%%%%%%%%%%
x = rand(N,1);
v = abs(rand(num_groups,N));

cvx_begin
    variable Z_cvx(num_groups,N)
    Mixed_norm = 0;
    for ii=1:num_groups
        Mixed_norm = Mixed_norm + norm(Z_cvx(ii,:),Inf);
    end
    minimize( rho*0.5*sum_square( vec(Z_cvx-G*diag(x)+v)) + lambda*Mixed_norm)
cvx_end

v=vec(v);

Z_cvx=Z_cvx(:);
[ z_admm ] = zupdate_1inf_overlap( x,1,G,1,1,v,rho,params, 1 );


max(abs(z_admm-Z_cvx))