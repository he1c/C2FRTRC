function Y = SNN_L1( T, Mask, I, option)
% dataname could be a string or actual data
% lambda1 = lambdaS * r * rRatio
% mode = RPCA mode (see rpca_for_tensor)
% mu1 = mu1fac * std(T(:))
addpath utils\SNN_L1\lightspeed
Omega=find(Mask(:)==1);
b=T(Omega);
N = ndims(T);
mu1fac = option.mulfac; 
rRatio = option.rRatio;

% params.lambda = 0.06;
r = 1/sqrt( max(size(T)) ); %0.015; %0.09:-0.005:0.02;
mu = mu1fac*std(b);  %mu1fac*std(T(:));
flag_lambdaS = 1;
lambda = flag_lambdaS*r*rRatio;       %params.lambdaS*r/4;
max_iter = option.maxitr;
debug = option.debug;
stopc = option.stopc; %1e-5 for syn data analysis plots, 1e-3 for others

%%%%%%%%%% for PROPACK %%%%%%%%%%%%
% declare global var 'sv'
global sv;
global tmode;
global use_propack;
global curr_mu;
sv =  ceil(min(size(T)) * 0.1) * ones( 1, N );
use_propack = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Solve HO-RPCA with missing data
%       min_{X,E} \sum_i ||X_(i)||_* + \lambda1*||E||_1
%       s.t. A_c(X + E) = T_c
% 
% Reformulate as
%       min \sum_i ||Xi,(i)||_* + \lambda1*||E||_1
%       s.t.    Xi = Y
%               A_c(Y+E) = T_c
%
% data.b = observations (vector)
% params
% X, V are cell arrays of tensors.
%
% Algorithm: ADAL

X0 = tenzeros( size(T) );
N = ndims(T);
V0 = cell( 1, N );
for i = 1:N
    V0{i} = tenzeros( size(T) );
end

tic;
bnorm = norm(b);
N = length( size(X0) );
X = cell( 1, N );
U = cell( 1, N );
V = cell( 1, N+1 );
for i = 1:N
    X{i} = X0;
    V{i} = V0{i};
end
V{N+1} = zeros( size(b) );
Y = X0;

% lambdaS = [ 1/5, 1/5, 1/5, 1 ];
if flag_lambdaS == -1
    lambdaS = [ 1, 1, 0.1, 0.1 ];
%     lambdaS = [ 0.01, 0.01, 1 ];
    lambda = -lambda;
else
    lambdaS = ones( 1, N );
end
rel_err = 1;

% assign operators
Aprod = @(X_param)A_select( X_param, Omega );
Atprod = @(b_param)At_select( b_param, Omega, size(Y) );

for iter = 1:max_iter
    % solve X_i's
    for i = 1:N
        [X{i}, junk, U{i}] = tensor_shrinkage( Y+mu*V{i}, mu*lambdaS(i), i );
    end
    
    % solve E
    P = Atprod( b - Y(Omega) + mu*V{N+1} );
    E = shrinkage_t( P, mu*lambda );
    
    % solve Y
    Yprev = Y;
    tensum = ten_sum_all( X )/mu - ten_sum_all( V(1:N) );
    Y = tensum * mu / N;
    Y(Omega) = (tensum(Omega) + V{N+1} + (b-E(Omega))/mu) / (N/mu + 1/mu);
    
    % compute optimality stats
    pres = 0;
    tdiff = cell( 1, N+2 );
    for i = 1:N
        tdiff{i} = X{i} - Y;
        pres = pres + norm( tdiff{i} )^2;
    end
    tdiff{N+1} = Aprod(Y+E)-b;
    pres = pres + norm(tdiff{N+1})^2;     pres = sqrt(pres);
    Ynorm2 = norm(Y)^2;
    Ydiff = Y - Yprev;
    dres = N*Ydiff;  dres(Omega) = (N-1)*Ydiff(Omega);  dres = norm(dres);
    denom = N*Yprev;    denom(Omega) = (N-1)*Yprev(Omega);     denom = norm(denom);
    pres = pres / sqrt( N*Ynorm2 + bnorm );%norm(b-Aprod(Y))^2 );    % / ynorm;
    dres = dres / denom;
    
    rel_err_p = rel_err;
    err=double(Y)-T;
    rel_err = norm(err(:)) / norm(T(:));
    imp = abs(rel_err-rel_err_p) / rel_err_p;
       
    if max(pres, dres) < stopc || imp < stopc
        break;
    end
    
    % update Lagrange multipliers
    for i = 1:length(V)
        V{i} = V{i} - tdiff{i}/mu;
    end
    
    if debug
        err1=double(Y)-I;
        disp(norm(err1(:)))
    end
    
end
rmpath utils\SNN_L1\lightspeed
end

function AX = A_select( X, Omega )
% X could be a tensor
% AX is a vector

X = double(X);
AX = X(Omega);

end

function X = At_select( b, Omega, sz )

X = zeros( sz );
X(Omega) = b;
X = tensor(X);
end

