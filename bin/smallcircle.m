function h = smallcircle(x,y,r)
% hold on
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
color = rand(1,3);
h = plot(xunit, yunit,'Color',color);
% hold off