function G = true_solution(zeta,alpha,q)


zeta1 = zeta;
zeta2 = zeta1.*alpha';
zeta1 = zeta1.*alpha^-1;
%t = @(zeta) q^2*zeta;

%omega = @(a,b) (t(a)-b)*(t(b) - a)/(t(a)-a)/(t(b)-b);
%w1 = (zeta1-alpha)*omega(zeta1,alpha);
%w2 = -w1/(alpha*zeta1);
%alph = (alpha')^-1;
%w3 = (zeta1-alph)*omega(zeta1,alph);
%w4 = -w3/(zeta1*alph);

%G1 = -1/(4*pi)*log(abs(w1*w2/w3/w4));
%disp(G1);

N = 100;

p1 = P(zeta1,q,N);
p2 = P(zeta2,q,N);

G2 = -1/(2*pi)*(log(abs(alpha*p1./p2)));

G = G2;
end
