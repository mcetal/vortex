function res = P(zeta, q, N)
    
    R = ones(76,65);
    for k=1:N
    
            Q1 = q^(2*k);
        
            term1 = Q1.*zeta;
    
            term2 = Q1.*(zeta.^-1);
    
            term3  = 1 + q^(4*k);
    
            Q = term3 - term1 - term2;
            
            R = R.*Q;
    end
    
    res = R.*(1-zeta);
end
    