function [T_est,error] = tam(m,n,k,r,l)
%UNTITLED Summary of this function goes here
% 
%
% % parameter setting
% 
% the tensor is m * n * k  
%          m   = 25;       
%          n   = 25;
%          k = 20;
%   r   = 5;                           % the tubal-rank
%   l=15                               %the number of iterations
% low-tubal-rank tensor
%    T = rand(m,n,k);                  %a ranom tensor: m * n * k
%    T = t_svd_threshold(T,r);         %make it to be tubal-rank = r
%profile on
T = tprod(rand(m,r,k), rand(r,n,k));
%[m,n,k] = size(T);
T_f = fft(T, [], 3);
% observations
p=0.05:0.05:0.25;
error = zeros(1,5);
T_omega = zeros(m,n,k);

for ii = 1:1
    omega = rand(m,n,k) <= 0.5;
    T_omega = omega .* T;

    T_omega_f = fft(T_omega,[],3);
    omega_f = fft(omega, [], 3);
% X: m * r * k
% Y: r * n * k
%% Given Y, do LS to get X
    Y = rand(r, n, k);
    %Y= init(T_omega, m,r,k);
    
    %[U, Theta, V]=t_svd(T_omega);
    %Y = V(1:r, :, :);
    
    Y_f = fft(Y, [], 3);

% do the transpose for each frontal slice
    Y_f_trans = zeros(n,r,k);
    X_f = zeros(m,r,k);
    T_omega_f_trans = zeros(n,m,k);
    omega_f_trans = zeros(n,m,k);
for i = 1: k
     Y_f_trans(:,:,i) = Y_f(:,:,i)';
     T_omega_f_trans(:,:,i) = T_omega_f(:,:,i)';
     omega_f_trans(:,:,i) = omega_f(:,:,i)';
end

iter=1;
while iter <=l
    fprintf('Sampling--%f---Round--#%d\n', p(ii), iter);
    [X_f_trans] = alter_min_LS_one_step(T_omega_f_trans, omega_f_trans * 1/k, Y_f_trans);
    
    for i =1:k
        X_f(:,:,i) = X_f_trans(:,:,i)';
    end

    % Given X, do LS to get Y
    [Y_f] = alter_min_LS_one_step(T_omega_f, omega_f * 1/k, X_f);
    
    for i = 1: k
    Y_f_trans(:,:,i) = Y_f(:,:,i)';
    end
    
    iter = iter + 1;
end

% The relative error:
temp = 0;
X_est = ifft(X_f, [], 3); 
Y_est = ifft(Y_f, [], 3);
T_est = tprod(X_est, Y_est);

temp = T - T_est;   
error(ii) = norm(temp(:)) / norm(T(:));
end

%profile viewer

end

function [Y_f] = alter_min_LS_one_step(T_omega_f, omega_f, X_f)
% the target dimension: r * s * k

[m,n,k] = size(T_omega_f);
[m,r,k] = size(X_f);

Y_f = zeros(r, n, k);
%% prepare variables
%X_f_new = zeros(m*k, r*k);
%for i=1:k
%    X_f_new((i-1)*m + 1 : i*m, (i-1)* r + 1 : i*r ) = X_f(:,:,i);
%end

X_f_new = [];
for i=1:k
    X_f_new = blkdiag(X_f_new, X_f(:,:,i));  % make it block diagonal
end


%% we recover the lateral slice one-by-one.
residual = 0;
tensor_V = zeros(k * m, 1);
temp_Y_f = zeros(r * k, 1);
for i = 1:n
    for j = 1:k
        tensor_V((j-1)*m + 1 : j*m) = squeeze(T_omega_f(:,i,j)); 
    end
    
    omega_f_3D = zeros(k,k,m);
    omega_f_new = zeros(k*m, k*m);
    for j = 1:m
        % omega_f_3D(:,:,j) = circ(squeeze(omega_f(j,i,:)));
        temp = gallery('circul',reshape(omega_f(j,i,:),1,k));
        omega_f_3D(:,:,j) = temp.';
    end
    
    for a=1:k
        for b=1:k
            for c=1:m
                row = (a-1)*m + c;
                col = (b-1)*m + c;
                omega_f_new(row,col) = omega_f_3D(a,b,c);
            end
        end
    end
    
    %[temp_Y_f, restnorm] = lsqlin(sparse(omega_f_new) * sparse(X_f_new), tensor_V);
    temp = sparse(omega_f_new) * sparse(X_f_new);
    temp_Y_f = temp \ tensor_V;
    %residual = residual + sqrt(restnorm)
    for j = 1:k
        Y_f(:,i,j) = temp_Y_f((j-1)*r + 1 : j*r);
    end 
end
end
