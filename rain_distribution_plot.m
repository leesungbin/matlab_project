syms N(d)

r=20/6;  %mm /hr
mu=5;
N(d)=8000*d^mu*exp(-3.67/1.4*d);  %m^-3 mm^-1
int(N,d,0,8)
int(N*d,d,0,8)/int(N,d,0,8)

range=0:.1:10;
plot(range,N(range));
