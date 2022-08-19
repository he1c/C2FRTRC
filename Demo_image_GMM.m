clear

addpath(genpath(pwd))
rmpath utils\SNN_L1\lightspeed

I=imread('data\BSD\78019.jpg');
 
I=double(I(1:320,1:480,:))./255;
opt.I=I;
frame=1;
opt.framenum=frame;
option=[];

%% parameter settings
option.debug=0;
option.maxitr=100;
n_t=[4,4,4,5,4,4,5,6,3];

N_t=length(n_t);
d=ceil(N_t/2);
I_t=reshape(I,n_t);
n_I=size(I);

v1=0.001;
v2=0.25;
c=0.4;
p=0.3;
   
opt.p=p;
disp(p)

omega = find(rand(numel(I),1)<p);
Mask = zeros(size(I));
Mask(omega) = 1;
Mask_t=reshape(Mask,n_t);

time=[];

disp(c);

G=noisemix(numel(I),1,c,v1,v2,'gaussian');
I_n=I+reshape(G,size(I));
I_n(I_n>1)=1;
I_n(I_n<0)=0;
MissM=Mask.*I_n;
MissM_t=reshape(MissM,n_t);

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
disp(['TTNN_L1 : ' num2str(psnr(I_TTNN_L1,I)) ' , time : ' num2str(time(11))])

%% TRC-Lp
option=get_option_default(MissM,'TRNN_Lp',opt);
t1=tic;
Utr = Completion_TR_lp_eps(MissM, Mask, option, I);
I_TRNN_Lp = reshape(Ui2U(Utr), n_I);
time(12)=toc(t1);
disp(['TRNN_Lp : ' num2str(psnr(I_TRNN_Lp,I)) ' , time : ' num2str(time(12))])


