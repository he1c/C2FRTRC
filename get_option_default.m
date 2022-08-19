function option=get_option_default(MissM,alg,opt)

option=[];

option.I=opt.I;
framenum=opt.framenum;

if ndims(option.I)==3&&framenum~=1
    option.nChannel=1;
else
    option.nChannel=3;
end

option.stopc = 1e-3;
option.maxitr = 100;
option.debug = 0;

option.yita = 4; %huber2
option.sigmamin = 0.15; %huber0.05
option.estimator = 'cauchy';

switch alg
    
    case 'HQTRC'
%         n_t=size(MissM);
%         N_t=length(n_t);
%         d=ceil(N_t/2);
%         sk=[];
%         if framenum==1
%             q=0.2;
%         else
%             q=0.1;
%         end
%         for n=1:N_t
%             order=[n:N_t 1:n-1];
%             M=reshape(MissM,prod(n_t(order(1:d))),[]);
%             sk=[sk max(ceil((min(size(M)))*q*sqrt(opt.p)),floor(sqrt(framenum))*2)];
%         end
%        option.r=sk; 
        option.r=ceil(sqrt(opt.p*size(option.I,1)*size(option.I,2))/5)*ones(1,ndims(MissM));
        option.mu=1e-4;
        option.lambda=2;  
        option.alpha=1.1;
        option.maxitr = 50;
    
    case 'PHQTRC'
        if framenum>1
            option.dim=20;
        else
            option.dim=36;
        end
        if isfield(opt,'dim')
            option.dim=opt.dim;
        end
        option.search_window_size=2;     
        option.stride=floor(option.dim*0.8);
        K=(2*option.search_window_size+1)^2;
        patch_miss=zeros(option.dim,option.dim,option.nChannel,framenum,K);
        patch_miss=squeeze(patch_miss);
%         n_t=size(patch_miss);
%         N_t=length(n_t);
%         d=ceil(N_t/2);
%         sk=[];
%         q=0.1;
%         for n=1:N_t
%             order=[n:N_t 1:n-1];
%             M=reshape(patch_miss,prod(n_t(order(1:d))),[]);
%             sk=[sk max(ceil((min(size(M)))*q*sqrt(opt.p)),floor(sqrt(K))*2)];
%         end
%         option.r=sk; 
        option.r=floor(sqrt(opt.p*option.dim*option.dim)*framenum^(1/3)/2)*ones(1,ndims(patch_miss));
        option.mu=1e-4;
        option.lambda=2;  
        option.alpha=1.1;
        option.wsigma=0.3;
        option.w0=0.2;
        if framenum>1
            option.maxitr = 10;
        else
            option.maxitr = 20;
        end
        option.yita2=0.02;
        option.yita3=4;
        option.localonly=0;

              
    case 'TRNN'
        option.beta = 1e-5;
        option.alpha = 1e-4;
        
    case 'TNN'
        option.rho=0.01;
        option.alpha=1;
        option.size=size(MissM);
        option.myNorm= 'tSVD_1';
        
    case 'TRNN_L1'
        option.mu=1e-4;
        
    case 'TNN_L1'
        option.beta    = 0.05;
        option.tauttnn = 1.618;
        option.lambda  = 1.5/sqrt(max(size(MissM,1),size(MissM,2))*size(MissM,3));
        
    case 'TMAC'
        option.maxit = 300; 
        option.tol = 1e-5; % run to maxit by using negative tolerance
        option.Mtr = MissM; % pass the true tensor to calculate the fitting
        option.alpha_adj = 0;
        option.rank_adj = -1*ones(1,3);
        if framenum==1
            option.rank_min = 10*ones(1,3);
            option.rank_max = 10*ones(1,3);
        else
            option.rank_min = 20*ones(1,3);
            option.rank_max = 50*ones(1,3);
        end
        
    case 'SNN_L1'
        option.stopc  = 1e-4;
        option.mulfac = 100;
        option.rRatio = 5;
        
    case 'SNN_WST'
        option.stopc    = 1e-4;
        option.alpha    = 1;
        option.lambda   = 0.05*min([size(MissM,1) size(MissM,2)])/sqrt(max([size(MissM,1) size(MissM,2)]));
        option.lrchannel= option.nChannel;
        
    case 'TNTV'
        option.beta   = 0.1;
        option.gamma  = 0.01;
        option.lambda = ceil(framenum*0.7)/sqrt(opt.p*size(MissM,3)*max(size(MissM,1),size(MissM,2)));
        option.rho    = 400;
        option.tau    = 1.618;
        option.theta  = 50;
        option.xi     = 1.2;
        option.tol    = 5e-3; 
        option.itemax = 100;
        
    case 'TTNN_L1'
        option.Data_Size = size(MissM);
        option.max_tot = 10^-4;
        option.max_iter= 20;
        option.disp = 1;
        r=ceil(sqrt(sqrt(opt.p*size(option.I,1)*size(option.I,2))/5));
        if framenum==1
            option.r = [r 2];
        else
            option.r = [r r 2];
        end
        option.p = 1;
    
    case 'TRNN_Lp'
        option.Data_Size = size(MissM);
        option.max_tot = 10^-4;
        option.max_iter= 20;
        option.disp = 1;
        option.r=ceil(sqrt(sqrt(opt.p*size(option.I,1)*size(option.I,2))/5))*ones(1,ndims(MissM));
        option.p = 1;  


end

end