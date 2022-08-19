function patch_arr = get_patch_ind(idxTotalPatch, param, framenum)

[m, n, ~] = size(idxTotalPatch);
stride     =   param.stride;  
border_ext =  floor(param.border_ext/2); 

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
        patch_arr((j-1)*length(ind_ref_m) + i)  = ind_refpatch;
    end
end
end

