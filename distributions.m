hold off
% for mu=[0 1 2 5 10]
mu=2
    n0=8000;
    org=@(t) t.^(mu).*exp(-t) %mm/hr
    A=(4.1*0.6^-0.21); %cm-1
    Nt=n0*integral(org,0,100)./A^(mu+1)
    
%     nt=@(d) n0.*d.^(mu+1).*exp(-(3.67+mu)/1.4.*d)/((3.67+mu)/1.4)^(mu+1)
N = @(d) Nt.*d.^mu.*exp(-A.*d); %rain distribution

total_d=integral(N,0,8)
p=@(d) integral(N,0,d)/total_d; %rain distribution

range=[0:.01:10];
plot(range,(N(range)))
hold on
% end