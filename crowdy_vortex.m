close all
clear
theta = 0:2*pi/64:2*pi;
q = 0.25;
r = q:0.01:1;
zeta_1 = exp(1i*theta);
zeta_2 = q*exp(1i*theta);


th = repmat(zeta_1,size(r,2),1);
ro = repmat(r',1,size(zeta_1,2));
zeta_test = ro.*th;
%
% in conformally mapped plane
figure(1)
   plot(real(zeta_1),imag(zeta_1),'r')
   hold on
   plot(real(zeta_2),imag(zeta_2),'b')
   title('zeta')

disp('pick an alpha')
xy = ginput(1)
alpha = xy(1) + 1i*xy(2)
   plot(xy(1),xy(2),'r*')

z=inline('(zeta-alpha)./(abs(alpha)*(zeta-1/alpha))','zeta','alpha')
z1 = z(zeta_1,alpha);
z2 = z(zeta_2,alpha);
x1 = inline('(z+conj(z))./(1+abs(z).^2)','z')
x2 = inline('(z-conj(z))./(1i*(1+abs(z).^2))','z')
x3 = inline('(-1+abs(z).^2)./(1+abs(z).^2)','z')

%zet = @(z) (alpha - (abs(alpha).*z)./alpha)./(1 - abs(alpha).*z);
%
% in stereographic plane
%Load points 
%load xzeta_grid.dat
 %a = zeros(NX,NY);
 %a(:) = xzeta_grid(:);
 %xz_grid = a;
%load yzeta_grid.dat
 %a = zeros(NX,NY);
 %a(:) = yzeta_grid(:);
 %yz_grid = a;
    
 %xzeta_grid = zet(xz_grid);
 %yzeta_grid = zet(yz_grid);
 
 %zeta_grid = xzeta_grid + 1i.*yzeta_grid;
 xzeta_test = real(zeta_test);
 yzeta_test = imag(zeta_test);
 gtest = true_solution(zeta_test,alpha,q);
 xz_test = real(z(zeta_test, alpha));
 yz_test = imag(z(zeta_test,alpha));
 
figure(2)
   plot(real(z1),imag(z1),'r')
   hold on
   plot(real(z2),imag(z2),'b')
   title('z')
    vc = -.5:.1:5;
    contour(xz_test,yz_test,gtest,vc)
    hold on
   
  
% on sphere
figure(3)
   sphere
   colormap(gray)
   shading flat
   hold on
   %alpha(0.5)
   plot3(x1(z1),x2(z1),x3(z1),'r')
   hold on
   plot3(x1(z2),x2(z2),x3(z2),'b')
   title('sphere')
