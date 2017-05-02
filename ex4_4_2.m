[y1,y2]=meshgrid(-10:.5:10,-10:.5:10);
dy1=y2;
dy2=-5*y1-2*y2;
quiver(y1,y2,dy1,dy2)
hold on

for theta=[0:16]*pi/8
    init1=8*[cos(theta);sin(theta)];
    [t,y]=ode45(@dydx1,[0 20],init1);
    plot(y(:,1),y(:,2),'k')
    hold on
    
%     init2=1*[theta;-6];
%     [t,y]=ode45(@dydx1,[0 20],init2);
%     plot(y(:,1),y(:,2),'k')
%     hold on
end
axis([-6 6 -6 6])
hold off