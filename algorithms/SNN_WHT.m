function M=SNN_WHT(T,Mask,I,option)

sigmamin=option.sigmamin;
stopc=option.stopc;
qtmin=option.qtmin;
debug=option.debug;
itr=option.maxitr;
alpha=option.alpha;
channel=option.lrchannel;
yita=option.yita;
r=option.rank;

[n1,n2,n3]=size(T);

X=cell(3,1);
Y=cell(3,1);
permute_order=cell(3,1);

permute_order{1}=[1 2 3];
permute_order{2}=[2 1 3];
permute_order{3}=[3 2 1];

permute_order_new{1}=[1 2 3];
permute_order_new{2}=[2 1 3];
permute_order_new{3}=[3 2 1];

for i=1:1:channel   
    X{i}=zeros(n1,n2,n3);
end

M=zeros(n1,n2,n3);

flag=0;

CC=0;

nn=[n1 n2 n3];

normT=norm(T(:));

for k=1:1:itr
           
    J=Mask.*(T-M);
    AAA=J(:);
    AAA(AAA==0)=[];
    
    sigma=max(max(abs(quantile(AAA,qtmin)),abs(quantile(AAA,1-qtmin)))*yita,sigmamin);


    for i=1:1:channel
        
        MaskMCC=Mask.*exp(-(M-T).^2/sigma^2);
        
        Diff_J=MaskMCC.*(M-T);
    
        Y{i}=X{i}-1/alpha*Diff_J;
              
        C = permute(Y{i},permute_order{i});
        
        Y_unfold = reshape(C,[],size(C,2),1);
        
        [U,S,V]=svd(Y_unfold,'econ');
        
        diagS_new=diag(S);
        
        diagS_new(r+1:end)=0;
        
        X_unfold=U*diag(diagS_new)*V';
        
        X_temp=reshape(X_unfold,nn(permute_order{i}(1)),nn(permute_order{i}(2)),nn(permute_order{i}(3)));
        
        X{i}= permute(X_temp,permute_order_new{i});
        
        M=0;
    
        for j=1:1:channel   
           M=M+X{j};       
        end
          
    end

    MaskMCCsq=sqrt(MaskMCC);
    
    CC_pre=CC; 
         
    CC=MaskMCCsq.*(T-M); 

    err=norm(CC(:)-CC_pre(:))/normT;
    
    if abs(err)<stopc
        break;
    end
    
    if debug
        err1=I-M;    
        disp(norm(err1(:)))
    end
        
end

end
