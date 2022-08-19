function Y=TRunfold(X,k,L)

order=[k:3 1:k-1];

Xp=permute(X,order);

if L==1
    Y=reshape(Xp,size(Xp,1),size(Xp,2)*size(Xp,3));
elseif L==2
    Y=reshape(Xp,size(Xp,1)*size(Xp,2),size(Xp,3));
elseif L==3
    Y=reshape(Xp,size(Xp,1)*size(Xp,2)*size(Xp,3),1);
end

end