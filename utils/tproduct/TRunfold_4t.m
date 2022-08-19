function Y=TRunfold_4t(X,k,L)

order=[k:4 1:k-1];

Xp=permute(X,order);

if L==1
    Y=reshape(Xp,size(Xp,1),size(Xp,2)*size(Xp,3)*size(Xp,4));
elseif L==2
    Y=reshape(Xp,size(Xp,1)*size(Xp,2),size(Xp,3)*size(Xp,4));
elseif L==3
    Y=reshape(Xp,size(Xp,1)*size(Xp,2)*size(Xp,3),1);
end

end