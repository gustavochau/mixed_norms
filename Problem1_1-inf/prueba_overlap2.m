clc;
clear;
close all;

N = 30;
num_groups = 50;
y = rand(N,1);

% create indices
rng(8)
for ii=1:num_groups
    p = randperm(N);
    g_index{ii} = p(1:randi(N/2,1,1));
end

lambda = 0.1;

cvx_begin
    variable x_cvx(N,1)
    Mixed_norm = 0;
    for ii=1:num_groups
        Mixed_norm = Mixed_norm + norm(x_cvx(g_index{ii}),Inf);
    end
    minimize( 0.5*sum_square( y - x_cvx) + lambda* Mixed_norm )
cvx_end

%%%%%%%%%%%%%%%%%
%%% pseudo admm
%%%%%%%%%%%%%%%%%%

rho=1;



G = zeros(num_groups,N);
for ii=1:num_groups
    G(ii,g_index{ii})=1;
end

H = [];
for ii=1:num_groups
%     G_aux = [G_aux; spdiag(G(ii,:))];
    H = [H; spdiags(G(ii,:)', 0, N, N)];
end    
% 
% for jj=1:50
% 
%     x = cvx_x( y,rho,z,G,v , N);
%     disp(['x update:' num2str(0.5*sum_square(vec(y-x))+lambda*compute_mixed_norm(z,1,inf)+0.5*rho*sum_square(vec(G*x-z)))])
%     %% z update
%     z_anterior = z;
%     z = cvx_z( y,rho,x,G,v,num_groups,N , lambda);
%     disp(['z update:' num2str(0.5*sum_square(vec(y-x))+lambda*compute_mixed_norm(z,1,inf)+0.5*rho*sum_square(vec(G*x-z)))])    
%     v = v+G*x-z;
%     disp(['v update:' num2str(0.5*sum_square(vec(y-x))+lambda*compute_mixed_norm(z,1,inf)+0.5*rho*sum_square(vec(G*x-z)))])    
% 
%     
%     costo(jj) = 0.5*sum_square(vec(y-x))+lambda*compute_mixed_norm(z,1,inf)+0.5*rho*sum_square(vec(G*x-z));
% 
%     
%     primal_res(jj) = norm(G*x - z,'fro');
%     dual_res(jj) = norm(rho*G'*(z-z_anterior),'fro');
%     
% end
% x_admm=full(diag(x));
% max(abs(x_admm-x_cvx))


z = zeros(num_groups*N,1);
v = zeros(num_groups*N,1);


costo_func = @(x1,z1) 0.5*sum_square(y-x1)+lambda*compute_mixed_norm( reshape(z1,N,num_groups)',1,inf)+0.5*rho*sum_square(H*x1-z1);

for jj=1:50

    x = cvx_x_2( y,rho,z,H,v,N );
    disp(['x update:' num2str(costo_func(x,z))])
    %% z update
    z_anterior = z;
    z = cvx_z_2( rho,x,G,v,num_groups,N ,lambda);
    disp(['z update:' num2str(costo_func(x,z))])
    v = v+H*x-z;
    disp(['v update:' num2str(costo_func(x,z))])

    
    costo(jj) = costo_func(x,z);

    
    primal_res(jj) = norm(H*x - z,'fro');
    dual_res(jj) = norm(rho*H'*(z-z_anterior),'fro');
    
end
x_admm=x;
max(abs(x_admm-x_cvx))

