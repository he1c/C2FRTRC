function [L,E]=TNN_L1(MissM,Mask,option,I)

beta=option.beta;
tau=option.tauttnn;
lambda=option.lambda;
maxitr=option.maxitr;
debug=option.debug;
stopc=option.stopc;

[n1,n2,n3]=size(MissM);

%% initial
Mask_C=~Mask;
L=zeros(n1,n2,n3);
E=L;
Z=L;

%% compute
for k=1:1:maxitr
    
    L_pre=L;

    M=MissM+Mask_C.*(L+E-1/beta*Z);

    [L,~,~]=prox_tnn(M+1/beta*Z-E,1/beta);

    M=MissM+Mask_C.*(L+E-1/beta*Z);

    H=M+1/beta*Z-L;

    E=sign(H).*max(abs(H)-lambda/beta,0);

    Z=Z-tau*beta*(L+E-M);
    
    J=L-L_pre;
    err=norm(J(:))/norm(MissM(:));
    
    if abs(err)<stopc
        break;
    end

    if debug
        err1=L-I;        
        nerr=norm(err1(:));
        disp(10*log10(n1*n2*n3/norm(nerr(:))^2))
    end
    
    
end