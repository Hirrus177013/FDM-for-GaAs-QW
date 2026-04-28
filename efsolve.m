function Ef = efsolve(Nd, A, beta, E, Ef_old)
    % EFSOLVE: Halley's Method (Newton-Raphson bậc 3)
    % Tự động tương thích với mảng E (15x1)
    
    Ef = Ef_old;
    tol = 1e-15;       
    max_iter = 50;     
    
    for iter = 1:max_iter
        exp_term = exp(beta * (Ef - E));
        
        % Chặn tràn bộ nhớ (Overflow)
        exp_term(isinf(exp_term)) = 1e300; 
        
        f_val = A * sum(log(1 + exp_term)) - Nd;
        f_prime = A * beta * sum(exp_term ./ (1 + exp_term));
        f_double_prime = A * (beta^2) * sum(exp_term ./ ((1 + exp_term).^2));
        
        numerator = 2 * f_val * f_prime;
        denominator = 2 * (f_prime^2) - f_val * f_double_prime;
        
        dEf = - (numerator / denominator);
        
        % Giới hạn bước nhảy tối đa là 0.5 eV
        max_step = 0.5 * 1.602e-19; 
        if abs(dEf) > max_step
            dEf = sign(dEf) * max_step;
        end
        
        Ef = Ef + dEf;
        
        if abs(dEf) < tol || abs(f_val) < 1e-10 * Nd
            break;
        end
    end
end