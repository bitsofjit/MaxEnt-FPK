clear all
clc
close all
x = linspace(-1.0,1.0,100);
y = x;
h = x(2)-x(1);
eps = 5;
alpha = eps*h;
x0 = 0.25; 
y0 = 0.3;
[X,Y] = meshgrid(x,y);
for i=1:length(x)
    for j=1:length(y)
        r(i,j) = sqrt((x(i)-x0)^2 + (y(j)-y0)^2);
    end
end
figure
rad_bas = exp((-alpha^2/h^2)*r.^2);
surfc(X,Y,rad_bas)
Ra = 0.5;
eps0 = 1.0e-1;
beta = -log(eps0)/Ra^2;
figure
rad_bas2 = exp(-beta*r.^2);
surfc(X,Y,rad_bas2)



