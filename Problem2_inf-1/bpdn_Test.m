close all;
clear;
clc;


Ig = double( imread('Imagenes/peppers.png') ) / 255;
[Nrows, Ncols]= size(Ig);

%%%%%%%%%%%%%%%%%%%%%%%%
%     Normalize        %
%%%%%%%%%%%%%%%%%%%%%%%%
Normalize = @(x) (x - min(x(:)))/(max(x(:)) - min(x(:)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             noisy images            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigma = 0.05;

% -------------
rng(5);
In = imnoise(Ig,'gaussian', 0, (sigma^2)*max(Ig(:)) );




%-----------------------------------------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          BPDN           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 0.5 || Phi u - b ||_2^2 + \lambda || u ||_1


% Regularization parameter (recommended values)
% 
lambda = 0.7;

% DCT Dictionary
dctLevels = 4;
K = @(x) dct2D_phi(x(:), Nrows, Ncols, dctLevels);
KT = @(x) dct2D_phiT(x(:), Nrows, Ncols, dctLevels);
Phi = {K, KT};
disp('===============')
disp('myADMM')
disp('===============')
opts.maxiter = 10; % máximo número de iteraciones
opts.rho0 = 1; % rho inicial
opts.tol = [1E-4 1E-4]; % [tolerancia absoluta  tolerancia_relativa]
opts.parrho = [5 1.5 1.5];
opts.u0 =randn(dctLevels*dctLevels*Nrows*Ncols,1); %solución inicial
opts.rhoopt = 'fix'; % opción de rho
% [u_admmF, statsF] = myADMM(Phi, In(:), lambda, opts); 
% IrecF = reshape(K(u_admmF(:,end)), Nrows, Ncols );
% figure, imshowpair(In,IrecF,'montage')

% 
% %%%%%%%% general ADMM
num_ele = dctLevels*dctLevels*Nrows*Ncols;
opts.x0 =opts.u0;
opts.z0 =opts.u0;

params.K=K;
params.KT=KT;
params.In = In(:);
params.lambda = lambda;
% solvers{1} = @(x,z,A,B,c,v,rho,params) bpdn_xupdate(x,z,A,B,c,v,rho,params);
% solvers{2} = @(x,z,A,B,c,v,rho,params) bpdn_zupdate(x,z,A,B,c,v,rho,params);
solvers{1} = @bpdn_xupdate;
solvers{2} = @bpdn_zupdate;
matrices{1} = speye(num_ele); %A
matrices{2} = -speye(num_ele); %B
matrices{3} = zeros(num_ele,1); %c
matrices{4} = matrices{1}; %A^T
func_costo = @(x,z) 0.5 * (norm(K(x)-In(:),2))^2 + lambda*(norm(x,1));
disp('===============')
disp('general ADMM')
disp('===============')
[u_iter,~] = generalADMM(matrices, func_costo, solvers,  params, opts);
Irec_gen = reshape(K(u_iter(end).x), Nrows, Ncols );
figure, imshowpair(In,Irec_gen,'montage')

% for ii=1:opts.maxiter
%     norm(u_admmF(:,ii)-u_iter(ii).x)
%     norm(statsF.z_iter(:,ii)-u_iter(ii).z)
% end