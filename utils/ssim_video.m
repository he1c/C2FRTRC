function y=ssim_video(I,I_recover)

if ndims(I)>3

    n_I=size(I);

    y=0;

    for i=1:1:length(n_I(end))
        y=y+ssim(I(:,:,:,i),I_recover(:,:,:,i));
    end

    y=y/length(n_I(end));

else
    
    y=ssim(I,I_recover);
end

end