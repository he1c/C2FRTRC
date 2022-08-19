
% clear; 
% close all; 
% 
% addpath('prox_operators');
% addpath('dataset');
% addpath('quality_assess');
% addpath('tensor_toolbox');
% addpath('my_main');

%% data claire
% load('Cl')
% Ltr = Cl./max(Cl(:));
% Ltr=Ltr(:,:,1:10);

%%
[n1 n2 n3] = size(Ltr);
 sr = 0.6;  %% sampling ratio
 level = 0.1;  %% impulse noise ratio
 sigma = 0.000001;  %% variance of Gaussian noise

Bn = imnoise(Ltr,'salt & pepper',level);  %% add impulse noise
noisy = Bn + sigma*randn([n1,n2,n3]);  %% noisy tensor with full observations

%% sampling ratio and Omega

fprintf('Sampling ratio = %0.8e\n',sr);
temp = randperm(n1*n2*n3);
kks = round((sr)*n1*n2*n3);
Omega = zeros(n1,n2,n3); 
Omega(temp(1:kks)) = 1;  %% index set of observations

Y = noisy.*Omega;  %% observed tensor

% %% parameters 
% 
% beta = 0.1;
% gamma = 0.01;
% lambda = 7/sqrt(sr*n3*max(n1,n2));
% rho1 = 400;
% tau = 1.618;
% tol = 5e-3;
% itemax = 100;
% theta = 50;
% xi = 1.2;
% 
% opts.theta  = theta;
% opts.xi     = xi;
% opts.tau    = tau;
% opts.rho    = rho1;
% opts.lambda = lambda;
% opts.beta   = beta;
% opts.gamma  = gamma;
% 
% opts.tol    = tol;   
% opts.itemax = itemax;
% 
% 
% %% main loop FFT
% fprintf('FFT \n');
% tic;
% [Lf S k eta beta] = LowSparN(Y, Omega, opts);
% toc;
%  %% Printf PSNR and SSIM
%     Error = norm(Lf(:)-Ltr(:))/norm(Lf(:));
%     fprintf('Relative error = %0.8e\n',Error);
%     PSNR = psnr(Ltr(:),Lf(:));
%     fprintf('PSNR = %0.8e\n',PSNR);
%     SSIM = ssim3d(Lf*255,Ltr*255);
%     fprintf('SSIM = %0.8e\n',SSIM);
    
    
%% DCT
 fprintf('DCT \n');

%% paraters            
opts.beta   = 0.1;
opts.gamma  = 0.01;
opts.lambda = 3/sqrt(sr*n3*max(n1,n2));
opts.rho    = 400;
opts.tau    = 1.618;
opts.tol    = 5e-3; 
opts.itemax = 100;

Ld = LowSparNDCT(Y, Omega, opts);

opts.rho = 1000;
O = Unfoldtntv(Ld,[n1 n2 n3],3);
[U,~,~] = svd(O,'econ');
Lu = LowSparNU(U, Y, Omega, opts);


%% Printf PSNR and SSIM
Error = norm(Lu(:)-Ltr(:))/norm(Lu(:));
fprintf('Relative error = %0.8e\n',Error);
PSNR = psnr(Ltr(:),Lu(:));
fprintf('PSNR = %0.8e\n',PSNR);
SSIM = ssim3d(Lu*255,Ltr*255);
fprintf('SSIM = %0.8e\n',SSIM);


