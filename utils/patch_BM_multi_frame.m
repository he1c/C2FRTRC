function [pos_arr, error_arr] = patch_BM_multi_frame(extractPatch, extractPatchMask, extractPatchref, extractPatchrefMask, idxTotalPatch, param, framenum)

[m, n, ~] = size(idxTotalPatch);
search_window_size    =   param.search_window_size;
K=param.K;
stride     =   param.stride;  
border_ext =  param.border_ext; 
% row / col index of search window left-up pos
% ind_ref_m  =  1:stride:m; 
% ind_ref_n   =  1:stride:n;

ind_ref_m =  border_ext:stride:m-border_ext; 
if ind_ref_m(end)<m-border_ext
    ind_ref_m=[ind_ref_m m-border_ext];
end
ind_ref_n  =  border_ext:stride:n-border_ext;
if ind_ref_n(end)<n-border_ext
    ind_ref_n=[ind_ref_n n-border_ext];
end


for  i  =  1 : length(ind_ref_m)
    for  j  =  1 : length(ind_ref_n)
        %// row / col
        row    =   ind_ref_m(i);
        col     =   ind_ref_n(j);
        ind_refpatch = idxTotalPatch(row, col, framenum);
        %// search window range <--> all searchable patch indices
        ind_search_m = max(row-search_window_size, 1) : min(row+search_window_size, m);
        ind_search_n = max(col-search_window_size, 1) : min(col+search_window_size, n);
        idx   =   idxTotalPatch(ind_search_m, ind_search_n,:);
        idx   =   idx(:);
        %// central patch
        meanCandidateMiss  =   extractPatchref(:,ind_refpatch);
        meanCandidateMissMask = extractPatchrefMask(:,ind_refpatch);
        ind_ob=meanCandidateMissMask~=0;
        %// all the patches in the region
        patchCandidateMiss = extractPatch(ind_ob,idx);
        patchCandidateMissMask = extractPatchMask(ind_ob,idx);
        %// distance: Euclidean & sorting
        meanCandidateMiss=meanCandidateMiss(ind_ob,:);
        meanCandidateMissMask=meanCandidateMissMask(ind_ob,:);
        disMiss = single(patchCandidateMiss) - single(meanCandidateMiss);
        disMissMask = patchCandidateMissMask.*(repmat(meanCandidateMissMask,1,length(idx)));
        dis = disMiss.*single(disMissMask);
        metric=   sum(dis.^2)./sum(single(disMissMask));
        [error, ind] =   sort(metric);
        %// take tensorSize-largest
        pos_arr(:, (j-1)*length(ind_ref_m) + i)    = idx(ind(1:K));
        error_arr(:, (j-1)*length(ind_ref_m) + i)  =  error(1:K);
    end
end
end

