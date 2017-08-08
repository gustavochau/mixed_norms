clc;
clear;
close all;

N = 6;
num_groups = 4;
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


z = rand(num_groups,N)*10;
v = 0*(rand(num_groups,N));

G = zeros(num_groups,N);
for ii=1:num_groups
    G(ii,g_index{ii})=1;
end

G_aux = [];
for ii=1:num_groups
%     G_aux = [G_aux; spdiag(G(ii,:))];
    G_aux = [G_aux; spdiags(G(ii,:)', 0, N, N)];
end    

for jj=1:5

    x = cvx_x( y,rho,z,G,v , N);
    disp(['x update:' num2str(0.5*sum_square(y-x)+lambda*compute_mixed_norm(z,1,inf)+0.5*rho*sum_square(vec(G*diag(x)-z)))])
    %% z update
    z_anterior = z;
    z = cvx_z( y,rho,x,G,v,num_groups,N , lambda);
    disp(['z update:' num2str(0.5*sum_square(y-x)+lambda*compute_mixed_norm(z,1,inf)+0.5*rho*sum_square(vec(G*diag(x)-z)))])    
    v = v+G*diag(x)-z;
    disp(['v update:' num2str(0.5*sum_square(y-x)+lambda*compute_mixed_norm(z,1,inf)+0.5*rho*sum_square(vec(G*diag(x)-z)))])    

    
    costo(jj) = 0.5*sum_square(y-x)+lambda*compute_mixed_norm(z,1,inf)+0.5*rho*sum_square(vec(G*diag(x)-z));

    
    primal_res = norm(G*diag(x) - z)
    dual_res = norm(rho*G'*(z-z_anterior))
    
end
x_admm=x;
max(abs(x_admm-x_cvx))

