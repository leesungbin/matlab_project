r0=1.06/1000;
r=[];
theta=[];
z=[];
count=0;
tc=0;
k=1;
while k<100
for i=[1:8000]
    r(i)=rand()*100*r0;
    theta(i)=rand()*2*pi;
    %z(i)=rand()*v*dt;
    if r(i)<2*r
        count=count+1;
    end
end
tc(k)=count;
count=0;
k=k+1;
end
