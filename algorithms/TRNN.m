function X=TRNN(T,Mask,I,option)

Ts=size(T);
N=ndims(T);

beta=option.beta;
alpha=option.alpha*ones(1,N);
stopc=option.stopc;
maxitr=option.maxitr;
debug=option.debug;

d=ceil(N/2);

J=size(T);

X=T;

for k=1:1:N
    M{k}=T;
    Y{k}=zeros(Ts);
end

for kk=1:1:maxitr
    
    X_pre=X;
    
    %% update M
    for k=1:1:N
        A=TRunfold_t(X,k,d)+1/beta*TRunfold_t(Y{k},k,d);
        [U,S,V]=svd(A,'econ');
        Sv=diag(S);
        Sv=max(Sv-alpha(k)/beta,0);
        DA=U*diag(Sv)*V';
        M{k}=TRfold_t(DA,k,J,N);
    end
    
    %% update X
    X=0;
    for k=1:1:N
        X=X+(1-Mask).*(M{k}-1\beta*Y{k})/N;
    end
    X=X+Mask.*T;
    
    %% update Y
    for k=1:1:N
        Y{k}=Y{k}+beta*(X-M{k});
    end
    
    E=X-X_pre;
    err=norm(E(:))/norm(T(:));
    
    if kk>5&&err<stopc
        break;
    end
    
    if debug&&mod(kk,10)==0
        err1=X-I;
        fprintf('Iter %.0f, Diff %.2f\n',kk,norm(err1(:)));
    end
end

end