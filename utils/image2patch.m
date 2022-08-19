function X = image2patch(image, dim)
%MODULE_IM2PATCH Summary of this function goes here
% Goal : decompose the 3D tensor (e.g., video) into 2D patches
[m, n, k] = size(image);
XX = cell(dim^2,k);
kk  = 0;
for i  = 1:dim
    for j  = 1:dim
        kk   =  kk+1;
        for p = 1:k      
            blk  =  image(i : m-dim+i, j : n-dim+j, p);
            XX{kk,p}=blk(:)';
        end
    end
end

XX=reshape(XX,[],1);
X=cell2mat(XX);

end


