function Utr = Completion_TR_norm(U_Omega, P_Omega, para)
    % U_Omega: I1 * I2 *... *In tensor with missing entries filled by 0
    % P_Omega: I1 * I2 *... *In binary observation tensor
    % para:    para.max_iter, para.max_tot, para.r
    %
    % Utr:     Decomposition term
    
    % Initialization
    Utr = TR_Initialization0(U_Omega,  para.r);
    dummy   = Utr{end};
    
    d = length(para.r);    
    tot     =1;
    iter    =0;
    
    while(tot>=para.max_tot && iter <=para.max_iter)
        %tic
        for i=1:d
            Utr{i} = Updata_U1_temp(TensPermute(P_Omega, i), Utr([i:d, 1:i-1]),TensPermute(U_Omega, i));
            if i<d
                U_temp=reshape(TensPermute(Utr{i},2),size(Utr{i},2),[]);
                for j=1:1:size(U_temp,2)
                    U_temp(:,j)=U_temp(:,j)./norm(U_temp(:,j));
                end
                    Utr{i}=TensPermute(reshape(U_temp,size(Utr{i},2),size(Utr{i},1),size(Utr{i},3)),3);
            end
        end
        %toc
        % Change of the last term as error
        tot = norm(T2V(Utr{d} - dummy))/norm(T2V(dummy));
        dummy = Utr{end};
        iter =iter+1;
        if para.disp ==1
            disp(['At iteration ', num2str(iter), ', the last term change is ', num2str(tot)]);
        end
    end
end