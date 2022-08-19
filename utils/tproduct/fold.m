function A=fold(X,n3)

[n1,n2]=size(X);

n1_n=floor(n1/n3);

A=zeros(n1_n,n2,n3);

for i=0:1:n3-1
    A(:,:,i+1)=X(n1_n*i+1:n1_n*(i+1),:);
end

end