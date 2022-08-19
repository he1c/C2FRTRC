function Y=TRfold_4t(X,n1,n2,n3,n4,k)

if k==1
    Y=reshape(X,[n1,n2,n3,n4]);
elseif k==2
    Y=reshape(X,[n2,n3,n4,n1]);
    Y=permute(Y,[4,1,2,3]);
elseif k==3
    Y=reshape(X,[n3,n4,n1,n2]);
    Y=permute(Y,[3,4,1,2]);
elseif k==4
    Y=reshape(X,[n4,n1,n2,n3]);
    Y=permute(Y,[2,3,4,1]);
end

end