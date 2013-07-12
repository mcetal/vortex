close all; clear;
NX = 1001; NY = 1000;
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
 load xxi_grid.dat
 a = zeros(NX,NY);
 a(:) = xxi_grid(:);
 xxi_grid = a;
 load yxi_grid.dat
 a = zeros(NX,NY);
 a(:) = yxi_grid(:);
 yxi_grid = a;
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
contour(xzeta_grid,yzeta_grid,ugrid,50)
hold on
geo_stereo
figure(4)
contour(xxi_grid, yxi_grid, ugrid, 100)
hold on
%figure(5)
%plot(xxi_grid, yxi_grid)
%hold on
%geo_stereo
%axis ([-15 7 -15 15])