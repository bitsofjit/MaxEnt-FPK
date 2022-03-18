clc;
clear all;
% close all;

format long g;
filenamen = 'solntn.dat';
fileIDn = fopen(filenamen,'r');

% Reading all the lines:
formatSpec = '%f %f %f';
soln = fscanf(fileIDn,formatSpec,[3,inf]);
soln = soln';

X = soln(:,1);
Y = soln(:,2);
Un = soln(:,3);

tri = delaunay(X',Y');

% How many triangles are there?
[r,c] = size(tri);
disp('# of traiangles resulting from Delaunay: ')
disp(r)

% Plot it with TRISURF
figure(1)
h = trisurf(tri, X', Y', Un');
axis vis3d
axis normal
% axis([0 inf 0 inf 0 inf]);
l = light('Position',[-50 -15 29]);
set(gca,'CameraPosition',[10 -1 10]);
view (50,30);
lighting phong
shading interp
colorbar EastOutside
camlight
grid on
light;
title('@ time t = t_n');

filename0 = 'solnt0.dat';
fileID0 = fopen(filename0,'r');
sol0 = fscanf(fileID0,formatSpec,[3,inf]);
sol0 = sol0';

X = sol0(:,1);
Y = sol0(:,2);
U0 = sol0(:,3);

tri = delaunay(X',Y');
% How many triangles are there?
[r,c] = size(tri);
disp('# of traiangles resulting from Delaunay: ')
disp(r)

% Plot it with TRISURF
figure(2)
title('At time t = 0');
h0 = trisurf(tri, X', Y', U0');
axis vis3d
axis normal
l = light('Position',[-50 -15 29]);
set(gca,'CameraPosition',[10 -1 10]);
view (50,30);
lighting phong
shading interp
colorbar EastOutside
camlight
grid on
light; 
title('@ time t = 0');


