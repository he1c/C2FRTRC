function A=fold_cell(X,r_n,n3)

[n1,n2]=size(X);

r=cumsum(r_n);
r=[0 r];

A=cell(1,n3);

for i=0:1:n3-1
    A{i+1}=X(r(i+1)+1:r(i+2),:);
end

end