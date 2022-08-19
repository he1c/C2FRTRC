function A=unfold_cell(X)

[n1,n2,n3]=size(X);

A=cell(1,n3);

for i=1:1:n3
    A{i}=squeeze(X(:,:,i));
end

end