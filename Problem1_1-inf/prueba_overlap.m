clc;
clear;
close all;

N = 30;
num_groups = 400;
y = rand(N,1);

% create indices
% rng(5)
for ii=1:num_groups
    p = randperm(N);
    g_index{ii} = p(1:randi(N/2,1,1));
end

lambda = 0.5;

cvx_begin
    variable x_cvx(N,1)
    Mixed_norm = 0;
    for ii=1:num_groups
        Mixed_norm = Mixed_norm + norm(x_cvx(g_index{ii}),Inf);
    end
    minimize( 0.5*sum_square( y - x_cvx) + lambda* Mixed_norm )
cvx_end


G = zeros(num_groups,N);
for ii=1:num_groups
    G(ii,g_index{ii})=1;
end

cvx_begin
    variable x_g(N,1)
    Mixed_norm2 = 0;
    Gx = G*diag(x_g);
    for ii=1:num_groups
        Mixed_norm2 = Mixed_norm2 + norm(Gx(ii,:),Inf);
    end
    minimize( 0.5*sum_square( y - x_g) + lambda* Mixed_norm2 )
cvx_end

max(abs(x_g-x_cvx))


% % 
% [x_admm,stats,u_iter] = admm_1_inf_overlap(y,g_index,lambda);
% 
% max(abs(x_admm-x_cvx))

