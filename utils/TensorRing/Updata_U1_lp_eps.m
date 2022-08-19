function U1 = Updata_U1_lp_eps(Omega, U, X_Omega, p)
    % Input
    % Omega:    Observe idx  I1 * ... * In
    % X_Omega:  Missing data I1 * ... * In
    % U:        {U1,..., Un}, each Ui is a 3 mode tensor
    %
    % Output
    % U1:       update U1 to solve:
    % U1 =      argmin ||Omega \circ f(yU2..Ud) - X_Omega||_F^2

    B = TCP(U(2:end));              % r1 * (I2...In) * rn
    
    n = length(U);
    P = tens2mat(Omega,   1, 2:n);  % I_1 * (I_2...I_n)
    X = tens2mat(X_Omega, 1, 2:n);  % I_1 * (I_2...I_n)
    
    [rn, I1, r1] = size(U{1});
    U1 = nan([rn, I1, r1]);
    
    for i=1:I1
        A = tens2mat(TensPermute(B(:, (P(i,:)==1), :),2), 1, [2,3]);
        y = (X(i, (P(i,:)==1)))';
        U1(:, i, :) = reshape(ADMM_eps(A, y, 20, p), [rn, 1, r1]);
    end

end 