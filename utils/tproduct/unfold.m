function A=unfold(X)

[n1,n2,n3]=size(X);

A=[];

for i=1:1:n3
    A=[A;squeeze(X(:,:,i))];
end

end