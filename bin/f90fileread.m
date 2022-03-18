clc;
clear all;
close all;
filename = 'matlabin.dat';

fileID = fopen(filename,'r');

% Reading first few three lines:
numnod= fscanf(fileID,'%d',1);
numq2 = fscanf(fileID,'%d',1);
Ra = fscanf(fileID,'%f',1);
% Ra = 1.5*Ra;

formatSpec = '%d %f %f';
sizenod = [3 numnod];
nodcoord = fscanf(fileID,formatSpec,sizenod);
nodcoord = nodcoord';

sizequad = [3 numq2];
quadcoord = fscanf(fileID,formatSpec,sizequad);
quadcoord = quadcoord';

fclose(fileID);

disp('number of nodes: ')
disp(numnod)
disp('number of quadrature points: ')
disp(numq2)
disp('nodal coordinate array extent: ')
disp(size(nodcoord))
disp('quadrature points coordinate array extent: ')
disp(size(quadcoord));
disp('Support Radius: ')
disp(Ra)

% Determining and plotting support nodes for a specific Gauss point:

norm = zeros(numq2,numnod);

for j=1:numq2
    for i=1:numnod
        norm(j,i) = sqrt((nodcoord(i,2)-quadcoord(j,2))^2 + (nodcoord(i,3)-quadcoord(j,3))^2);
    end
end

nghbr = []; v = [];

for j=1:numq2
    counter = 0;  
    for i=1:numnod
        if (norm(j,i)-Ra) < 100*eps
            counter = counter + 1;
            v(j,counter) = i;
        end
    end
    nghbr = [nghbr; counter];
end

% These following 2 lines will make the plot messy. SO not using:
strValues1 = strtrim((num2str([nodcoord(:,1)],'(%d)')));
strValues2 = strtrim((num2str([quadcoord(:,1)],'(%d)')));

% The Gauss point no at non-convergence:
gp = input('Specify the Gauss point # where non-convergence occurs: ');

gp_nghbr_no = nghbr(gp);
[row,column,gp_support_nodes] = find(v(gp,:));

% Writing the Gauss point number on plot:
strValues3 = strtrim((num2str([quadcoord(gp,1)],'(%d)')));

% The support nodes for the above Gauss Points:
strValues4 = strtrim((num2str([nodcoord(gp_support_nodes,1)],'(%d)')));

% Plotting:
figure
% Don't use the following line- makes the plot messy
% text(nodcoord(:,2),nodcoord(:,3),strValues1,'VerticalAlignment','bottom');
plot(nodcoord(:,2),nodcoord(:,3),'ko','MarkerFaceColor','k','Markersize',6);
hold on;
text(nodcoord(gp_support_nodes,2),nodcoord(gp_support_nodes,3), ...
    strValues4,'VerticalAlignment','bottom','Fontsize',10);

hold on;
plot(quadcoord(:,2),quadcoord(:,3),'r.','MarkerFaceColor','r','Markersize',1);
% Don't use the following line- makes the plot messy
% text(quadcoord(:,2),quadcoord(:,3),strValues2,'VerticalAlignment','bottom');

% Just pointing out the un-cool Gauss point separately:
plot(quadcoord(gp,2),quadcoord(gp,3),'ro','MarkerFaceColor','r','Markersize',6);

axis('equal')

for i=gp_support_nodes
    h = smallcircle(nodcoord(i,2),nodcoord(i,3),Ra);
end
hold off

figure
plot(nodcoord(:,2),nodcoord(:,3),'ko','MarkerFaceColor','k','Markersize',6);
hold on;
plot(quadcoord(:,2),quadcoord(:,3),'r.','MarkerFaceColor','r','Markersize',6);