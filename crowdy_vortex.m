close all
clear
theta = 0:2*pi/64:2*pi;
q = 0.25;
zeta_1 = exp(1i*theta);
zeta_2 = q*exp(1i*theta);
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
%
% in stereographic plane
figure(2)
   plot(real(z1),imag(z1),'r')
   hold on
   plot(real(z2),imag(z2),'b')
   title('z')
%
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
