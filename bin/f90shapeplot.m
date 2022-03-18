clc;
clear all;
% close all;
format long g;
name = 'phi.dat';
fid = fopen(name,'r');
formatSpec = '%f %f %f';
sizeofarray = [3 inf];
matrixdata = fscanf(fid,formatSpec,sizeofarray);
matrixdata = matrixdata';

XX = matrixdata(:,1);
YY = matrixdata(:,2);
PHI = matrixdata(:,3);

% figure
% plot3(XX,YY,PHI,'r*')
% axis([-1 1 -1 1 0 1])
% grid on;
% title('MinXEnt Basis Functions');
% xlabel('X')
% ylabel('Y') 
% zlabel('MinXEnt Basis Function (Point plots)')
% hold off

% Little triangles
% The solution is to use Delaunay triangulation. Let's look at some
% info about the "tri" variable.
tri = delaunay(XX',YY');
figure
plot(XX',YY','r.','MarkerSize', 1)
hold on

% How many triangles are there?
[r,c] = size(tri);
disp('# of traiangles resulting from Delaunay: ')
disp(r)

% Plot it with TRISURF
h = trisurf(tri, XX', YY', PHI');
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
title('MinXEnt Basis Functions (Surface Plot)');
xlabel('X')
ylabel('Y') 
zlabel('MinXEnt Basis Function')
