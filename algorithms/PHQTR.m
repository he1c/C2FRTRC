function Xhat=PHQTR(MissM,Mask,Mhat,option)

%% parameter setting
dim=option.dim;
search_window_size=option.search_window_size;
K=(2*search_window_size+1)^2;
nChannel=option.nChannel;
stride=option.stride;
I=option.I;
wsigma_min=option.wsigma;
w0_max=option.w0;
yita2=option.yita2;
yita3=option.yita3;
localonly=option.localonly;


%% pre-processing
[m, n, ~]=size(MissM);
if nChannel==3
    framenum=size(MissM,4);
elseif nChannel==1
    framenum=size(MissM,3);
end

MissM = padarray(MissM,[search_window_size search_window_size],'symmetric','both');
Mask = padarray(Mask,[search_window_size search_window_size],'symmetric','both');
Mhat = padarray(Mhat,[search_window_size search_window_size],'symmetric','both');
I = padarray(I,[search_window_size search_window_size],'symmetric','both');
numpatch=prod([dim,dim,nChannel,framenum,K]);

MissM(Mask==0)=-1;
   
%% initial
m=m+2*search_window_size;
n=n+2*search_window_size;
m2=m-dim+1;
n2=n-dim+1;
total_patch_size = [m2, n2];

ind_ref_m =  search_window_size+1:stride:m2-search_window_size; 
if ind_ref_m(end)<m2-search_window_size
    ind_ref_m=[ind_ref_m m2-search_window_size];
end
ind_ref_n = search_window_size+1:stride:n2-search_window_size;
if ind_ref_n(end)<n2-search_window_size
    ind_ref_n=[ind_ref_n n2-search_window_size];
end

ind_patch=[];
for i = 1 : length(ind_ref_m)
    for j = 1 : length(ind_ref_n)
        row = ind_ref_m(i);
        col = ind_ref_n(j);
        ind_patch = [ind_patch sub2ind(total_patch_size, row, col)];
    end
end

%% patch jitter
blk_arr_n=zeros(K-1,length(ind_patch));
for i=1:1:length(ind_patch)
    blk_arr_n(:,i)=patch_jitter(total_patch_size,ind_patch(i),search_window_size);
end
blk_arr_n=[ind_patch;blk_arr_n];

patch_hat=cell(1,size(blk_arr_n,2));

for kk=1:size(blk_arr_n,2)

    arr=blk_arr_n(:,kk);
    patch_miss=zeros(dim,dim,nChannel,framenum,K);  
    patch_ref=zeros(dim,dim,nChannel,framenum,K);  
    for i=1:1:K
        [xx,yy]=ind2sub(total_patch_size,arr(i));
        xxn=xx:min(xx+dim-1,m);
        yyn=yy:min(yy+dim-1,n);
        patch_miss(1:length(xxn),1:length(yyn),:,:,i)=MissM(xxn,yyn,:,:); 
        patch_ref(1:length(xxn),1:length(yyn),:,:,i)=Mhat(xxn,yyn,:,:); 
    end
    patch_miss=squeeze(patch_miss);
    patch_ref=squeeze(patch_ref);
    
    mask_patch=double(patch_miss~=-1);
    patch_comb=patch_miss.*mask_patch+patch_ref.*(1-mask_patch);
    errpatch=(patch_comb-patch_ref).*mask_patch;
    J=errpatch(:);
    J=J(mask_patch(:)==1);
    p=length(J)/numpatch;
    wsigma=max(yita2/max(abs(quantile(J,0.75)),abs(quantile(J,0.25))),wsigma_min);
    W=exp(-errpatch.^2/2/wsigma^2).*mask_patch;
    w0=min(yita3*max(abs(quantile(J,0.75)),abs(quantile(J,0.25))),w0_max);
    r=floor(sqrt(p*dim*dim)*framenum^(1/3)/2)*ones(1,ndims(patch_miss));

    if localonly==1
        JJ=patch_miss(:);
        JJ=JJ(mask_patch(:)==1);
        Xinitial=mean(JJ);
        patch_hat{kk}=RWTRR(patch_miss,mask_patch,mask_patch,option,Xinitial,r);
    else
        W=w0*(1-mask_patch)+W.*mask_patch;
        Xinitial=mean(patch_ref(:));
        patch_hat{kk}=RWTRR(patch_comb,double(patch_comb~=-1),W,option,Xinitial,r);
    end

end

%% aggreagate
if framenum==1

    I_hat=zeros(m,n,nChannel);
    C_hat=zeros(m,n,1);

    for kk=1:1:size(blk_arr_n,2)
          for ii=1:1:K
            [xx,yy]=ind2sub(total_patch_size,blk_arr_n(ii,kk));
            xxn=xx:min(xx+dim-1,m);
            yyn=yy:min(yy+dim-1,n);
            I_hat(xxn,yyn,:)=I_hat(xxn,yyn,:)+patch_hat{kk}(1:length(xxn),1:length(yyn),:,ii);   
            C_hat(xxn,yyn)=C_hat(xxn,yyn)+ones(length(xxn),length(yyn));  
          end

    end

    
elseif nChannel>1&&framenum>1
    
    I_hat=zeros(m,n,nChannel,framenum);
    C_hat=zeros(m,n,1);

    for kk=1:1:size(blk_arr_n,2)
          for ii=1:1:K
            [xx,yy]=ind2sub(total_patch_size,blk_arr_n(ii,kk));
            xxn=xx:min(xx+dim-1,m);
            yyn=yy:min(yy+dim-1,n);
            I_hat(xxn,yyn,:,:)=I_hat(xxn,yyn,:,:)+patch_hat{kk}(1:length(xxn),1:length(yyn),:,:,ii);   
            C_hat(xxn,yyn)=C_hat(xxn,yyn)+ones(length(xxn),length(yyn));  
          end

    end
    
elseif nChannel==1&&framenum>1
    
    I_hat=zeros(m,n,framenum);
    C_hat=zeros(m,n,1);

    for kk=1:1:size(blk_arr_n,2)
          for ii=1:1:K
            [xx,yy]=ind2sub(total_patch_size,blk_arr_n(ii,kk));
            xxn=xx:min(xx+dim-1,m);
            yyn=yy:min(yy+dim-1,n);
            I_hat(xxn,yyn,:,:)=I_hat(xxn,yyn,:,:)+patch_hat{kk}(1:length(xxn),1:length(yyn),:,ii);   
            C_hat(xxn,yyn)=C_hat(xxn,yyn)+ones(length(xxn),length(yyn));  
          end

    end

       
end

Xhat=I_hat(search_window_size+1:end-search_window_size,search_window_size+1:end-search_window_size,:,:)./C_hat(search_window_size+1:end-search_window_size,search_window_size+1:end-search_window_size);

Xhat(isnan(Xhat))=0; 

Xhat=squeeze(Xhat);

end
