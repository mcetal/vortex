function res = P(zeta, q, N)
    k = 1:N;
    
    Q1 = q.^(2*k);
    
    Q2 =  1 - zeta*Q1;
    
    Q = Q2.*(1 - zeta^-1*Q1);
    
    disp(Q);
    
    res = (1-zeta)*Q;
end
    