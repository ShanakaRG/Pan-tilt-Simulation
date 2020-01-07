clear all;
close all;
clc;
crc = @(r1,r2,p,a) [r1*cosd(a + p);  r2*sind(a + p)];
t = linspace(0, 360, 361); % Time
init_pos = 0;              % Initial Position (°)
r1 = 4;
r2 = 2;
locus = crc(r1,r2, init_pos, t)

figure(1)
for k1 = 2:length(t)
    x= locus(1,k1)
    y=locus(2,k1)
    plot(x, y , 'r*')
    axis([-5  5    -5  5])
    grid
    axis square
    drawnow    

end


%  x1=5;
%  x2=-5;
%  y1=0;
%  y2=0;
%  e=0.8
%  a = 1/2*sqrt((x2-x1)^2+(y2-y1)^2);
%  b = a*sqrt(1-e^2);
%  t = linspace(0,2*pi);
%  X = a*cos(t);
%  Y = b*sin(t);
%  w = atan2(y2-y1,x2-x1);
%  x = (x1+x2)/2 + X*cos(w) - Y*sin(w);
%  y = (y1+y2)/2 + X*sin(w) + Y*cos(w);
%  plot(x,y,'r-')
%  axis equal