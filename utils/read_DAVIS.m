clear

video_name='flamingo';
SamplePath1 =  ['data\DAVIS\' video_name '\'];  
fileExt = '*.jpg';  

files = dir(fullfile(SamplePath1,fileExt)); 

frame =80;
%遍历路径下每一幅图像
I=[];
for i=1:1:frame
   fileName = strcat(SamplePath1,files(i).name); 
   image = imread(fileName);
   I(:,:,:,i) = image;  
end

I=imresize(I,0.25);

save(['data\DAVIS\' video_name '.mat'],'I')