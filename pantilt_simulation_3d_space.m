clear all;
close all;
clc;

a=1
%  while (2>1)
for i= pi/18:pi/180:14*pi/12
% Your two points
P1 = [0,0,0];

rho=400;
phi=cos(i)
theta=i+180
x=rho*sin(phi).*cos(theta);
y=rho*sin(phi).*sin(theta);
z=rho*cos(phi);
P2 = [x,y,z]; 
pts = [P1; P2];
line(pts(:,1), pts(:,2), pts(:,3));
f=plot3(pts(:,1), pts(:,2), pts(:,3));
% f=plot(x)
drawnow
hold on;
grid on
box on


f = getframe(gcf);
[im,map] = rgb2ind(f.cdata,256,'nodither');
im(:,:,1,a) = rgb2ind(f.cdata,map,'nodither');
% colormap(f.colormap);
% k=image(f.cdata);
% imshow(k)

%frames(a) = getframe(f);
%imwrite(k,['Image' int2str(a), '.jpg']);
a=a+1 
%imwrite(im,map,'testAnimated.gif','gif','LoopCount',Inf,'DelayTime',1);
end

% 
% for i=0:pi/180:pi
% % Your two points
% P1 = [0,0,0];
% 
% rho=400;
% phi=sin(i)
% theta=i
% x=rho*sin(phi).*cos(theta)
% y=rho*sin(phi).*sin(theta)
% z=rho*cos(phi)
% P2 = [x,y,z]; 
% pts = [P1; P2];
% line(pts(:,1), pts(:,2), pts(:,3))
% plot3(pts(:,1), pts(:,2), pts(:,3))
% drawnow
% hold on;
% grid on
% box on
% end
% end
% h = animatedline;
% axis([0,4*pi,-1,1])
% 
% x = linspace(0,4*pi,1000);
% y = sin(x);
% for k = 1:length(x)
%     addpoints(h,x(k),y(k));
%     drawnow
% end


% theta=linspace(0,2*pi,40);
% phi=linspace(0,pi,40);
% [theta,phi]=meshgrid(theta,phi);
% rho=1;
% x=rho*sin(phi).*cos(theta);
% y=rho*sin(phi).*sin(theta);
% z=rho*cos(phi);
% mesh(x,y,z)
% axis equal
% grid on
% box on
% view([130,30])
% xlabel('x-axis')
% ylabel('y-axis')
% title('The graph of x^2+y^2+z^2=1.')


% syms u v
% x = cos(u)*sin(v);
% y = sin(u)*sin(v);
% z = cos(v)*sin(v);
% figure
% fsurf(x,y,z)
% axis equal
% syms t
% 
% Rx = [1 0 0; 0 cos(t) -sin(t); 0 sin(t) cos(t)]
% Ry = [cos(t) 0 sin(t); 0 1 0; -sin(t) 0 cos(t)]
% Rz = [cos(t) -sin(t) 0; sin(t) cos(t) 0; 0 0 1]
% xyzRx = Rx*[x;y;z];
% Rx45 = subs(xyzRx, t, pi/4);
% figure
% fsurf(Rx45(1), Rx45(2), Rx45(3))
% title('Rotating by \pi/4 about x, counterclockwise')
% axis equal


% v1=[0.2,0.3 4 ]
% v2=[-0.3,0.3,0.1]
% v=[v2;v1]
% plot3(v(:,1),v(:,2),v(:,3),'r')
% 
% grid on
% box on
% 
% 
% 
% R=10;
% Phi=linspace(-pi/12,pi/12);
% Theta=linspace(0,pi/4);
% [Phi,Theta]=meshgrid(Phi,Theta);
% [X,Y,Z]=sph2cart(Theta,Phi,R);
% % Z=R*sin(Phi);
% % X=R*cos(Phi).*cos(Theta);
% % Y=R*cos(Phi).*sin(Theta);
% figure
% mesh(X,Y,Z);
% grid on
% box on

% r = linspace(0,2*pi) ;
% th = linspace(0,2*pi) ;
% [R,T] = meshgrid(r,th) ;
% X = R.*cos(T) ;
% Y = R.*sin(T) ;
% Z = R ;
% surf(X,Y,Z)
% 
% 
% figure(); hold on;
% for theta = linspace(0, pi, 100) % Not exactly sure how you want to vary theta
%     [T, R] = meshgrid(linspace(0, theta, 100), 1:100);
%     [X, Y] = pol2cart(T,R); 
%     Z = 100 - R.^2; % Compute the surface of revolution
%     surf(X,Y,Z); % Plot the surface
%     pause(1); % Wait one second
% end

% for i=1:10
% A=[1 2 3+i];
% B=[1 2 4+i];
% 
% 
% 
% % theta=linspace(0,2*pi,40);
% % phi=linspace(0,pi,40);
% % [theta,phi]=meshgrid(theta,phi);
% % rho=1;
% % x=rho*sin(phi).*cos(theta);
% % y=rho*sin(phi).*sin(theta);
% % z=rho*cos(phi);
% 
% 
% plot3(A,B,'k-')
% hold on
% end