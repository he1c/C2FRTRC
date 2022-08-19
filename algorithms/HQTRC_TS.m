function X_t=HQTRC_TS(M_t,Mask_t,I,option)

Ts=size(M_t);
N=ndims(M_t);
No=N;
r=option.r;
stopc=option.stopc;
alpha=option.alpha;
maxitr=option.maxitr;
debug=option.debug;
estimator=option.estimator;
yita=option.yita;
sigmamin=option.sigmamin;

lambda=option.lambda;
mu=option.mu;

d=ceil(N/2);
%N=d; %only use first half of the unfolding matrices

J=size(M_t);

X_t=zeros(Ts);

for k=1:1:N
    G_t{k}=zeros(Ts);
end

X_t_pre=0;

err=0;

count=0;

XG_temp=cell(1,N);
Z_t=XG_temp;

for kk=1:1:maxitr
    
    %%==========update W_t==========
    E_t=Mask_t.*(M_t-X_t);
    E_v=E_t(:);
    E_v(E_v==0)=[]; 
    sigma=max(max(abs(quantile(E_v,0.75)),abs(quantile(E_v,0.25)))*yita,sigmamin);
    if strcmp(estimator,'correntropy')
        W_t=exp(-E_t.^2/sigma^2/2);
    elseif strcmp(estimator,'cauchy')
        W_t=sigma.^2/(sigma.^2+E_t.^2);
    elseif strcmp(estimator,'huber')
        T=abs(E_t)<=sigma;
        W_t=T.*1+~T.*(sigma./(abs(E_t)+1e-10));
    else 
        disp('Invailad M-estimator');
        X_t=[];
        break;
    end
    
    %tic
    
    %%==========update Z_tk==========
    for k=1:1:N
        order=[k:No 1:k-1];
        XG_temp{k}=permute(X_t-G_t{k}/mu,order);
        XG_temp{k}=reshape(XG_temp{k},prod(J(order(1:d))),[]);
%         [U,S,V]=svd(XG_temp{k},'econ');
%         Sn=diag(S);
%         Sn(r(k)+1:end)=0;
%         Z_temp=reshape(U*diag(Sn)*V',J(order));
%         Z_t{k}=ipermute(Z_temp,order);
    end
    
    
    [U,S,V]=cellfun(@(A)svd(A,'econ'),XG_temp, 'UniformOutput',false);
    Z_temp=cellfun(@(A,B,C)A(:,1:r(k))*B(1:r(k),1:r(k))*C(:,1:r(k))',U,S,V, 'UniformOutput',false);
    
    for k=1:1:N
      order=[k:No 1:k-1];
      Z_temp2=reshape(Z_temp{k},J(order));
      Z_t{k}=ipermute(Z_temp2,order);
    end
    
    %toc

    
    %==========update X_t==========
    L=0;
    for k=1:1:N
        L=L+Z_t{k}+1/mu*G_t{k};
    end
    L=L/N;
    
    X_t=L+lambda*W_t.*Mask_t./(lambda*W_t+1).*(M_t-L);
    
    
    
    %%==========update G_tk==========
    
    for k=1:1:N
        G_t{k}=G_t{k}+mu*(Z_t{k}-X_t);
    end
    
    E=X_t-X_t_pre;
    
    %norm(E(:))
    
    err_pre=err;
    
    err=norm(E(:))/norm(X_t_pre(:));
    
    if abs(err-err_pre)<stopc
        count=count+1;
        if count>2
            if debug
                err1=X_t-I;
                fprintf('Iter %.0f, Diff %.2f\n',kk,norm(err1(:)));
            end
            break;
        end
    else
        count=0;
    end
        
    X_t_pre=X_t;
    
    if debug&&mod(kk,1)==0
        err1=X_t-I;
        fprintf('Iter %.0f, Diff %.2f, ObjDiff %.5f\n',kk,norm(err1(:)),err-err_pre);
    end
   
    mu=mu*alpha;

end

end