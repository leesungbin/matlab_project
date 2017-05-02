h=2000;
rho_a=-0.0977*h+1.225;

r=1.1443/1000;
V=4/3*pi*r^3;
A=pi*r^2;
rho_w=1000;

Cd=0.47;
mass=V*rho_w;
g=9.80;
v=0;
% F=mass*g-1/2*Cd*rho_a*A*v^2;
F=0.1;

dt=0.01;
dh=0;

i=2;

while h>0
    v(i)=F(i-1)*dt/mass+v(i-1);
    dh(i)=(v(i-1)+v(i))/2*dt;
    h(i)=h(i-1)-dh(i);
    rho_a(i)=-0.0977*h(i)/1000+1.225;
    F(i)=mass*g-1/2*Cd*rho_a(i)*A*v(i)^2;%-rho_a(i)*V*g; 

    i=i+1;
end
t=0:dt:(length(v)-1)*dt;
plot(t,v,'k')

xlabel("t(s)")
ylabel("v(m/s)")

