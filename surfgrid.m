close all; clear;
NX = 500; NY = 501;
% nx = nphi, ny = nth+1
load igrid.dat
 a = zeros(NX,NY);
 a(:) = igrid(:);
 igrid = a;
load xgrid.dat
 a = zeros(NX,NY);
 a(:) = xgrid(:);
 xgrid = a;
load ygrid.dat
 a = zeros(NX,NY);
 a(:) = ygrid(:);
 ygrid = a;
load zgrid.dat
 a = zeros(NX,NY);
 a(:) = zgrid(:);
 zgrid = a;
load ugrid.dat
 a = zeros(NX,NY);
 a(:) = ugrid(:);
 ugrid = a;
load xzeta_grid.dat
 a = zeros(NX,NY);
 a(:) = xzeta_grid(:);
 xzeta_grid = a;
load yzeta_grid.dat
 a = zeros(NX,NY);
 a(:) = yzeta_grid(:);
 yzeta_grid = a;
load ugrid.dat
 a = zeros(NX,NY);
 a(:) = ugrid(:);
 ugrid = a;
%
%
figure(1)
surf(xgrid,ygrid,zgrid,igrid)
   a=colormap(gray);
   new(1:64,:)=a(64:-1:1,:);
   colormap(new);
   caxis([0 3])
   shading flat
figure(2)
surf(xgrid,ygrid,zgrid,ugrid)
   shading flat
   hold on
   geo_3d
   caxis ([-.5 5])
figure(3)
vc = [-.5:.1:5];
contour(xzeta_grid,yzeta_grid,ugrid,vc)
hold on
