function [u, MSE] = ADMM(A, x, mu, p)
e = randn(size(x));
Lambda = randn(size(x));
MSE = [];
maxiter = 80;
g_diff = @(e_i,y_i) e_i - y_i +p/mu*e_i^(p-1);
g_diff2 = @(e_i,y_i) e_i - y_i -p/mu*(-e_i)^(p-1);
g = @(e_i,y_i) 0.5*(e_i - y_i)^2 + 1/mu*abs(e_i)^p;
for iter = 1 : maxiter
    u = A\(e - Lambda/mu + x);
    y = A*u+ Lambda/mu - x;
    
    if p ==1
        y_sign = y;
        y_sign( y_sign<0 )=-1;
        y_sign( y_sign>0 )=1;
        y_hat = abs(y) - 1/mu*ones(size(y));
        y_hat(y_hat<0)=0;
        e = y_sign.*y_hat;
    elseif p<1
        beta_a = 2/mu*(1-p)^(1/(2-p));
        h_a = beta_a + p/mu*beta_a^(p-1);
        for i = 1 : size(y,1)
            if abs(y(i)) <= h_a
                e(i) = 0;
            else
                beta = beta_a;
                for j = 1 : 20
                    beta = abs(y(i)) -1/mu*p*beta^(p-1);
                end
                beta_star = beta;
                e(i) = sign(y(i))*beta_star;
            end
        end
    else
        for i = 1 : size(y,1)
            if y(i) >= 0
                a = 0;
                b = y(i);
                c = (a + b)/2;
                middle = g_diff(c,y(i));
                inner_iter = 1;
                while abs(middle) > 1e-8
                    if (middle > 0)
                        b = c;
                    else
                        a = c;
                    end
                    c = (a+b)/2;
                    middle = g_diff(c,y(i));
                    inner_iter = inner_iter + 1;
                    if inner_iter > 50
                        break
                    end
                end
                if g(c,y(i)) < g(0,y(i))
                    e(i) = c;
                else
                    e(i) = 0;
                end
            else
                a = y(i);
                b = 0;
                c = (a + b)/2;
                middle = g_diff2(c,y(i));
                inner_iter = 1;
                while abs(middle) > 1e-8
                    if (middle > 0)
                        b = c;
                    else
                        a = c;
                    end
                    c = (a+b)/2;
                    middle = g_diff2(c,y(i));
                    inner_iter = inner_iter + 1;
                    if inner_iter > 50
                        break
                    end
                end
                if g(c,y(i)) < g(0,y(i))
                    e(i) = c;
                else
                    e(i) = 0;
                end
            end
        end
    end
    
    Lambda = Lambda + mu*(A*u - e - x);
    
end


end