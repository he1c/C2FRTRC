function Y=TRfold(X,n1,n2,n3,k)

if k==1
    Y=reshape(X,[n1,n2,n3]);
elseif k==2
    Y=reshape(X,[n2,n3,n1]);
    Y=permute(Y,[3,1,2]);
elseif k==3
    Y=reshape(X,[n3,n1,n2]);
    Y=permute(Y,[2,3,1]);
end

end