%% 1.读取视频内容并显示
video_name='highway';

fid = fopen(['data\YUV\' video_name '_cif.yuv']); %读入YUV文件
row=144*2;col=176*2; %图像的高、宽
frames=200; %序列的帧数，当前只处理5帧图像
im_l = {};
im_l_ycbcr = zeros([row, col, 3, frames]);
im_ycbcr = zeros(row, col);

for frame=1:frames
 %读入文件 将yuv转换为rgb，并用imshow显示
    im_l_y = zeros(row,col); %Y
    for i1 = 1:row 
       im_l_y(i1,:) = fread(fid,col);  %读取数据到矩阵中 
    end

    im_l_cb = zeros(row/2,col/2); %cb
    for i2 = 1:row/2 
       im_l_cb(i2,:) = fread(fid,col/2);  
    end

    im_l_cr = zeros(row/2,col/2); %cr
    for i3 = 1:row/2 
       im_l_cr(i3,:) = fread(fid,col/2);  
    end

    %由于输入的yuv文件为4:2:0，所以CbCr要改变大小，
    %否则im_l_ycbcr(:, :, 2) =im_l_cb;会出现错误
    im_l_cb = imresize(im_l_cb, [row, col], 'bicubic');%改变图像的大小
    im_l_cr = imresize(im_l_cr, [row, col], 'bicubic');
       
    im_ycbcr(:, :, 1) = im_l_y;
    im_ycbcr(:, :, 2) = im_l_cb;
    im_ycbcr(:, :, 3) = im_l_cr;
    
    im_l_ycbcr(:,:,:,frame)=uint8(ycbcr2rgb(uint8(im_ycbcr)));
    
end

I=uint8(im_l_ycbcr);
save(['data\YUV\' video_name '.mat'],'I')