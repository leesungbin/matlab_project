r=1.1443/1000;
V=4/3*pi*r^3;
A=pi*r^2;
rho_w=1000;
rho_a=((-0.0977*2+1.225)+1.225)/2;
Cd=0.47;
mass=V*rho_w;
g=9.80;
v=0;
F=mass*g-1/2*Cd*rho_a*A*v^2;

dt=0.01;
h=2000;
dh=0;

i=2;
while h>0
    v(i)=F(i-1)*dt/mass+v(i-1);
    dh(i)=(v(i-1)+v(i))/2*dt;
    h(i)=h(i-1)-dh(i);
    F(i)=mass*g-1/2*Cd*rho_a*A*v(i)^2-rho_a*V*g;

    i=i+1;
end
t=0:dt:(length(v)-1)*dt;
plot(t,v,'k')

xlabel("t(s)")
ylabel("v(m/s)")


%integral Fdt = m*vf
syms y(t)
eqn=g-1/2*Cd*rho_a*A*y(t)^2/mass==diff(y,t); %-rho_a*V*g/mass
solution(t)=dsolve(eqn,y(0)==0)