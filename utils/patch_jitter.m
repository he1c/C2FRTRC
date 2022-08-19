function pos_arr = patch_jitter(totalpatchsize, patchind, search_window_size)

[row,col]=ind2sub(totalpatchsize,patchind);
step=1;
K=(2*search_window_size+1)^2;
row_new=row;
col_new=col;
xxn=[];
yyn=[];
while(length(xxn)<K)
    step=-step;
    xxn=[xxn' row_new:sign(step):row_new+step (row_new+step)*ones(1,abs(step))]';
    yyn=[yyn' col_new*ones(1,abs(step)) col_new:sign(step):col_new+step]';
    step=sign(step)*(abs(step)+1);
    row_new=xxn(end);
    col_new=yyn(end);
    xxn(end)=[];
    yyn(end)=[];
end

idx = sub2ind(totalpatchsize, xxn(1:K), yyn(1:K));
idx(idx==patchind)=[];
pos_arr = idx;

end


