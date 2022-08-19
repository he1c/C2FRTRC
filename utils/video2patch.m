function X = video2patch(video, dim)
%MODULE_IM2PATCH Summary of this function goes here
% Goal : decompose the 3D tensor (e.g., video) into 2D patches
if ndims(video)<2
    disp('video dim must be 2, 3 or 4');
    X=[];
    return;
end
[m, n, k, numFrame] = size(video);
%numFrame=size(video,ndims(video));
N   =   m-dim+1;
M   =   n-dim+1;
L     =   N*M;                            % total #patches
X     =   zeros(dim^2*k, L, numFrame, class(video));        % n
for idxFrame = 1 : numFrame
    kk  = 0;
    for i  = 1:dim
        for j  = 1:dim
           kk   =  kk+1;
           if k==1
                blk  =  video(i : m-dim+i, j : n-dim+j, idxFrame);
                X(kk, :, idxFrame) =  blk(:)';
            elseif k==3
                for p=1:k
                    blk  =  video(i : m-dim+i, j : n-dim+j, p, idxFrame);
                    X(kk+dim^2*(p-1), :, idxFrame) =  blk(:)';
                end
           end
        end
    end
end
end

