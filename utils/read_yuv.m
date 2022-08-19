%% 1.��ȡ��Ƶ���ݲ���ʾ
video_name='highway';

fid = fopen(['data\YUV\' video_name '_cif.yuv']); %����YUV�ļ�
row=144*2;col=176*2; %ͼ��ĸߡ���
frames=200; %���е�֡������ǰֻ����5֡ͼ��
im_l = {};
im_l_ycbcr = zeros([row, col, 3, frames]);
im_ycbcr = zeros(row, col);

for frame=1:frames
 %�����ļ� ��yuvת��Ϊrgb������imshow��ʾ
    im_l_y = zeros(row,col); %Y
    for i1 = 1:row 
       im_l_y(i1,:) = fread(fid,col);  %��ȡ���ݵ������� 
    end

    im_l_cb = zeros(row/2,col/2); %cb
    for i2 = 1:row/2 
       im_l_cb(i2,:) = fread(fid,col/2);  
    end

    im_l_cr = zeros(row/2,col/2); %cr
    for i3 = 1:row/2 
       im_l_cr(i3,:) = fread(fid,col/2);  
    end

    %���������yuv�ļ�Ϊ4:2:0������CbCrҪ�ı��С��
    %����im_l_ycbcr(:, :, 2) =im_l_cb;����ִ���
    im_l_cb = imresize(im_l_cb, [row, col], 'bicubic');%�ı�ͼ��Ĵ�С
    im_l_cr = imresize(im_l_cr, [row, col], 'bicubic');
       
    im_ycbcr(:, :, 1) = im_l_y;
    im_ycbcr(:, :, 2) = im_l_cb;
    im_ycbcr(:, :, 3) = im_l_cr;
    
    im_l_ycbcr(:,:,:,frame)=uint8(ycbcr2rgb(uint8(im_ycbcr)));
    
end

I=uint8(im_l_ycbcr);
save(['data\YUV\' video_name '.mat'],'I')