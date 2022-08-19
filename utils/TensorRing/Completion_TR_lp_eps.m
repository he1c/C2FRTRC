function Utr = Completion_TR_lp_eps(U_Omega, P_Omega, para, Data)
    % U_Omega: I1 * I2 *... *In tensor with missing entries filled by 0
    % P_Omega: I1 * I2 *... *In binary observation tensor
    % para:    para.max_iter, para.max_tot, para.r
    %
    % Utr:     Decomposition term
    
    % Initialization
    Utr = TR_Initialization(U_Omega,  para.r);
    
    d = length(para.r);    
    tot     =1;
    iter    =0;
    Data_Size = para.Data_Size;
    while(tot > para.max_tot && iter <=para.max_iter)
        for i=1:d
            Utr{i} = Updata_U1_lp_eps(TensPermute(P_Omega, i), Utr([i:d, 1:i-1]),TensPermute(U_Omega, i),para.p);
        end
        % Change of the last term as error
        tot = norm(T2V(reshape(Ui2U(Utr), Data_Size) - Data./255))/norm(T2V(Data./255));
        iter =iter+1;
        
        if para.disp ==1
            disp(['At iteration ', num2str(iter), ', the last term change is ', num2str(tot)]);
        end
    end
end