clear

addpath(genpath(pwd))
rmpath utils\SNN_L1\lightspeed

load(['data\tempete.mat']); 
 
I=double(I)./255;

frame=30;

I=imresize(I(:,:,:,1:frame),[240 432]);
opt.I=I;
opt.framenum=frame;

%% parameter settings

p=0.8;
opt.p=p;
option.debug=0;
option.maxitr=100;
n_I=size(I);

v1=0.01;
v2=0.25;
c=0;

time=[];    


Mask=imbinarize(imread('watermark4.bmp'));
Mask=Mask(50:end-80,40:end-40,:);
Mask=double(imresize(Mask,[n_I(1),n_I(2)]));
Mask=repmat(Mask,1,1,1,frame);

G = zeros(n_I);
for f=1:1:frame
    for i=1:1:n_I(3)
        for j=1:1:n_I(1)
            if rand(1)<c
                G(j,:,i,f)=noisemix(n_I(2),1,1,v1,v2,'gaussian');
            else
                G(j,:,i,f)=noisemix(n_I(2),1,0,v1,v2,'gaussian');
            end
        end
    end
end
I_n=I+G;
MissM=Mask.*I_n;


for i=1:1:10
    imwrite(MissM(:,:,:,i),['results\videobmp\0000' num2str(i-1) '.bmp'])
end

for i=11:1:size(I,4)
    imwrite(MissM(:,:,:,i),['results\videobmp\000' num2str(i-1) '.bmp'])
end

for i=1:1:10
    imwrite(1-Mask(:,:,:,i),['results\videobmpmask\0000' num2str(i-1) '.bmp'])
end

for i=11:1:size(I,4)
    imwrite(1-Mask(:,:,:,i),['results\videobmpmask\000' num2str(i-1) '.bmp'])
end



% for nn=[5 8 10 14 15 16]
%     
%     load(['results\Group5\' num2str(nn) '_p1.mat']);
%     
%     disp(['I_OPTRC SSIM: ' num2str(psnr(I_OPTRC,I)) ' time:' num2str(time(end))]);
% 
%     I_STTN=zeros(size(I_OPTRC));
% 
%     for i=1:1:frame_num
%         Im=imread(['D:\STTN\' num2str(nn) '\result_' num2str(i) '.bmp']);
%         I_STTN(:,:,:,i)=double(imresize(Im,[n_o(1) n_o(2)]))/255;
%     end
% 
%     I_STTN=(1-Mask).*I_STTN+Mask.*MissM;
% 
% 
%     disp(['I_STTN SSIM: ' num2str(psnr(I_STTN,I)) ' time:' num2str(time(end))]);
% 
% 
%     I_FUSE=zeros(size(I_OPTRC));
% 
%     for i=1:1:frame_num
%         Im=imread(['D:\FuseFormer\' num2str(nn) '\result_' num2str(i) '.bmp']);
%         I_FUSE(:,:,:,i)=double(imresize(Im,[n_o(1) n_o(2)]))/255;
%     end
% 
%     I_FUSE=(1-Mask).*I_FUSE+Mask.*MissM;
% 
% 
%     disp(['I_FUSE SSIM: ' num2str(psnr(I_FUSE,I)) ' time:' num2str(time(end))]);
%     
%     save(['results\Group5\' num2str(nn) '_p1.mat']);
% 
% end
    