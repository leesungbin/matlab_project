Nt=8000;
mu=0:2;

D0=2.1; %median drop size - diameter
A=(3.67+mu)/D0;

d=0:.1:10;
for i=mu
    ND=Nt*d.^i.*exp(-(3.67+i)/D0.*d);
    plot(d,(ND))
    hold on
end
r=[0.85 1.125 2.45 4.35 5.1]
n=[1360 2720 2320 1280 320]
plot(r,n)
hold off
