clc;
clear all;
% close all;
format long g;
name = 'dphi.dat';
fid = fopen(name,'r');
formatSpec = '%f %f %f %f';
sizeofarray = [4 inf];
matrixdata = fscanf(fid,formatSpec,sizeofarray);
matrixdata = matrixdata';

XX = matrixdata(:,1);
YY = matrixdata(:,2);
dPHIdx = matrixdata(:,3);
dPHIdy = matrixdata(:,4);

% Little triangles
% The solution is to use Delaunay triangulation. Let's look at some
% info about the "tri" variable.
tri = delaunay(XX',YY');
% X derivative:
figure
plot(XX',YY','r.', 'MarkerSize', 1)
hold on

% How many triangles are there?
[r,c] = size(tri);
disp('# of traiangles resulting from Delaunay: ')
disp(r)

% Plot it with TRISURF
h = trisurf(tri, XX', YY', dPHIdx');
axis vis3d

% Clean it up:
axis normal
% axis off
l = light('Position',[-50 -15 29]);
set(gca,'CameraPosition',[10 -1 10]);
view (50,30);
lighting phong
shading interp
colorbar EastOutside
camlight
grid on
light; 

% colorbar;
title('X Derivative of maxent Basis Functions (Surface Plot)');
xlabel('X')
ylabel('Y') 
zlabel('ddX(max-ent) Basis Function')
hold off;

% Y derivative
figure
plot(XX',YY','r.', 'MarkerSize', 1)
hold on


[r,c] = size(tri);
disp('# of traiangles resulting from Delaunay: ')
disp(r)

% Plot it with TRISURF
h = trisurf(tri, XX', YY', dPHIdy');
axis vis3d

% Clean it up:
axis normal
% axis off
l = light('Position',[-50 -15 29]);
set(gca,'CameraPosition',[10 -1 10]);
view (50,30);
lighting phong
shading interp
colorbar EastOutside
camlight
grid on
light; 

title('Y Derivative of maxent Basis Functions (Surface Plot)');
xlabel('X')
ylabel('Y') 
zlabel('ddY(max-ent) Basis Function')
hold off;
