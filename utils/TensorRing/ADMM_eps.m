function [u, MSE] = ADMM_eps(A, x, mu, p)
e = randn(size(x));
Lambda = randn(size(x));
MSE = [];
maxiter = 40;
tolerance = 1e-10;
for iter = 1 : maxiter
    u = A\(e - Lambda/mu + x);
    y = A*u+ Lambda/mu - x;
    ep = sqrt(1-p)*max(norm(e,'inf'),norm(y,'inf'))+1e-2;
      
    if p<=1
        for m=1:1:length(e)
            if abs(y(m))<tolerance*100
                my_fun = @(x) p/mu*x*(x^2+ep^2)^(p/2-1) + x - y(m);
                low = - abs(y(m));
                high = abs(y(m));
                e(m) = bisection(my_fun, low, high, tolerance); 
            else
                e(m) = 0;
            end
        end
            
    end
    
    Lambda = Lambda + mu*(A*u - e - x);
    
end

end