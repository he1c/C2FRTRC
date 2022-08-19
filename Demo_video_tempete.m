clear

addpath(genpath(pwd))
rmpath utils\SNN_L1\lightspeed

load('data\tempete.mat'); 
 
I=double(I)./255;

frame=30;

I=imresize(I(:,:,:,1:frame),[144 180]);
opt.I=I;
opt.framenum=frame;

%% parameter settings
option.debug=0;
option.maxitr=100;
n_t=[3,3,4,4,3,3,4,5,3,5,6];

I_t=reshape(I,n_t);
n_I=size(I);

v1=0.0001;
v2=0.25;

time=[];

for c=0.5

    Mask=imbinarize(imread('data\watermark_word.bmp'));
    Mask=Mask(50:end-80,40:end-40,:);
    Mask=double(imresize(Mask,[n_I(1),n_I(2)]));
    Mask=repmat(Mask,1,1,1,frame);
    Mask_t=reshape(Mask,n_t);

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
    MissM_t=reshape(MissM,n_t);
    
    opt.p=sum(Mask(:))/numel(Mask);
    MissM3=reshape(MissM,n_I(1),n_I(2),[]);
    Mask3=reshape(Mask,n_I(1),n_I(2),[]);

    %% HQTRC
    option=get_option_default(MissM_t,'HQTRC',opt);
    t1=tic;
    M=HQTRC_TS(MissM_t,Mask_t,I_t,option);
    I_HQTRC_TS=reshape(M,n_I);
    time(1)=toc(t1);
    disp(['HQTRC : ' num2str(psnr(I_HQTRC_TS,I)) ' , time : ' num2str(time(1))])

    %% C2FRTRC
    option=get_option_default(MissM,'PHQTRC',opt);
    t1=tic;
    fprintf('patch size: %d, overlap: %d',option.dim,option.dim-option.stride);
    I_C2FRTRC=PHQTR(MissM,Mask,I_HQTRC_TS,option);
    time(2)=toc(t1);
    disp(['C2FRTRC : ' num2str(psnr(I_C2FRTRC,I)) ' , time : ' num2str(time(2))])

    %% LPRTRC
    option=get_option_default(MissM,'PHQTRC',opt);
    t1=tic;
    option.localonly=1;
    I_PHQTRC=PHQTR(MissM,Mask,I_HQTRC_TS,option);
    time(3)=toc(t1);
    disp(['LPRTRC : ' num2str(psnr(I_PHQTRC,I)) ' , time : ' num2str(time(3))])

    %% TNN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    option=get_option_default(MissM3,'TNN',opt);
    t1=tic;
    A = diag(sparse(double(Mask3(:)))); 
    b = MissM3(:);
    I_hat = TNN(A,b,option);
    I_TNN = reshape(I_hat,size(I));
    time(4)=toc(t1);
    disp(['TNN: ' num2str(psnr(I_TNN,I)) ' time:' num2str(time(4))]);

    %% TRNN
    option=get_option_default(MissM_t,'TRNN',opt);
    t1=tic;
    option.stopc=1e-3;
    M=TRNN(MissM_t,Mask_t,I_t,option);
    I_TRNN=reshape(M,n_I);
    time(5)=toc(t1);
    disp(['TRNN : ' num2str(psnr(I_TRNN,I)) ' , time : ' num2str(time(5))])

    %% TRNN_L1
    option=get_option_default(MissM,'TRNN_L1',opt);
    t1=tic;
    M=TRNN_L1(MissM,Mask,option,I);
    I_TRNN_L1=reshape(M,n_I);
    time(6)=toc(t1);
    disp(['TRNN_L1 : ' num2str(psnr(I_TRNN_L1,I)) ' , time : ' num2str(time(6))])

    %% TNN_L1
    option=get_option_default(MissM3,'TNN_L1',opt);
    t1=tic;
    [M,E]=TNN_L1(MissM3,Mask3,option,I);
    I_TNN_L1=reshape(M,n_I);
    time(7)=toc(t1);
    disp(['TNN_L1 : ' num2str(psnr(I_TNN_L1,I)) ' , time : ' num2str(time(7))])

    %% SNN_L1
    t1=tic;
    option=get_option_default(MissM3,'SNN_L1',opt);
    M = SNN_L1(MissM3,Mask3,I,option);
    I_SNN_L1=reshape(M.data,n_I);
    time(8)=toc(t1);
    disp(['SNN_L1 : ' num2str(psnr(I_SNN_L1,I)) ' , time : ' num2str(time(8))])

    %% SNN_WST
    t1=tic;
    option=get_option_default(MissM3,'SNN_WST',opt);
    M= SNN_WST(MissM3,Mask3,I,option);
    I_SNN_WST=reshape(M,n_I);
    time(9)=toc(t1);
    disp(['SNN_WST : ' num2str(psnr(I_SNN_WST,I)) ' , time : ' num2str(time(9))])

    %% TNTV
    option=get_option_default(MissM3,'TNTV',opt);
    t1=tic;
    Ld = LowSparNDCT(MissM3, Mask3, option);
    option.rho = 1000;
    O = Unfoldtntv(Ld,size(MissM3),3);
    [U,~,~] = svd(O,'econ');
    Lu = LowSparNU(U, MissM3, Mask3, option);
    I_TNTV=reshape(Lu,n_I);
    time(10)=toc(t1);
    disp(['TNTV : ' num2str(psnr(I_TNTV,I)) ' , time : ' num2str(time(10))])

    %% TTC-L1
    option=get_option_default(MissM,'TTNN_L1',opt);
    t1=tic;
    Utr = Completion_TT(MissM, Mask, option, I);
    I_TTNN_L1 = reshape(Ui2U(Utr), n_I);
    time(11)=toc(t1);
    disp(['TTC_L1 : ' num2str(psnr(I_TTNN_L1,I)) ' , time : ' num2str(time(11))])

    %% TRC-Lp
    option=get_option_default(MissM,'TRNN_Lp',opt);
    t1=tic;
    Utr = Completion_TR_lp_eps(MissM, Mask, option, I);
    I_TRNN_Lp = reshape(Ui2U(Utr), n_I);
    time(12)=toc(t1);
    disp(['TRC-Lp : ' num2str(psnr(I_TRNN_Lp,I)) ' , time : ' num2str(time(12))])
    
end


