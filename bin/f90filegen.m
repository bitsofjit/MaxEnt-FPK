clc;
clear all;
close all
format long g;

% Mesh file generation:
nsdim = 2;
ndivl = 50;
ndivw = 50;
xx = linspace(-10.0,10.0,ndivl);
yy = linspace(-10.0,10.0,ndivw);
[X,Y] = meshgrid(xx,yy);
TRI = delaunay(X,Y);

triplot(TRI, X(:), Y(:));

data = [X(:),Y(:)];
data = data';

numnod = size(data,2);
ntriangle = size(TRI,1);

x = data(1,:);
y = data(2,:);

coord = zeros(3,numnod);
coord(1,:) = 1:1:numnod;
coord(2,:) = y;
coord(3,:) = x;

% Check:
figure
strValues1 = strtrim((num2str([coord(1,:)'],'(%d)')));
text(coord(2,:),coord(3,:),strValues1,'VerticalAlignment','bottom');
hold on;
plot(coord(2,:)',coord(3,:)','r.','Markersize',6);
hold off;

TRI = TRI';
ITER = 1000;
TOL = 1.E-014;
OBJ_FUNC = 'jaynes';
PRIOR = 'gaussian-rbf';
strdim = num2str(nsdim);
strn = num2str(numnod);

filename1 = [strdim 'd_' strn '_Nodes_Structured_mesh.dat'];

fileID1 = fopen(filename1,'w');
fprintf(fileID1,'%5s %5s\n','NSD','NODES');
fprintf(fileID1,'%d %d\n',nsdim, numnod);
fprintf(fileID1,'%10s\n','ELEMS');
fprintf(fileID1,'%d\n',ntriangle);
fprintf(fileID1,'%10s\n','COORD');
fprintf(fileID1,'%d %20.16f %20.16f\n',coord);
fprintf(fileID1,'%10s\n','CONNECTIVITY');
fprintf(fileID1,'%d %d %d\n',TRI);
fprintf(fileID1,'%10s %10s\n','ITER','TOL');
fprintf(fileID1,'%d %20.16f\n',ITER, TOL);
fprintf(fileID1,'%10s %10s\n','PRIOR','OBJ-FUNC');
fprintf(fileID1,'%10s %10s\n',PRIOR,OBJ_FUNC);
fclose(fileID1);

% Determining the boundary nodes in the domain [0, 1] X [0, 1]
bcnod = [];
for i=1:numnod
    if (coord(2,i) == min(x)) || (coord(2,i) == max(x))
        bcnod = [bcnod i];
    end
end

for i=1:numnod
    if (coord(3,i) == min(y)) || (coord(3,i) == max(y))
        bcnod = [bcnod i];
    end
end

bcnod = unique(bcnod);
lthu = length(bcnod);
% Check:
figure
strValues2 = strtrim((num2str([bcnod'],'(%d)')));
text(coord(2,bcnod)',coord(3,bcnod)',strValues2,'VerticalAlignment','bottom');
hold on;
plot(coord(2,bcnod), coord(3,bcnod), 'r.','Markersize',6);

for i=1:lthu
    for j=1:size(coord,2)
        if (coord(1,j) == bcnod(i))
            nodval(i) = 0.0; 
        end
    end
end
bcdata = [bcnod; nodval];

filename2 = [strdim 'd_' strn '_Nodes_Structured_ebc.dat'];

fileID2 = fopen(filename2,'w');
fprintf(fileID2,'%10s\n','NUMBER');
fprintf(fileID2,'%d\n',lthu);
fprintf(fileID2,'%10s %10s\n','NOD', 'VAL');
fprintf(fileID2,'%d %20.16f\n',bcdata);
fclose(fileID2);






