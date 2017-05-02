hold off

% N = @(d) 8000.*d.^5.*exp(-3.67./1.4.*d); %rain gamma distribution
N = @(d) 8000.*d.^5.*exp(-3.67./1.4.*d); %rain distribution
total_d=integral(N,0,8);
p=@(d) integral(N,0,d)/total_d; %rain distribution

drop_count=0;
dt=0.1;
g=9.80;
h0=2000;
Cd=0.47;
r_max=2.5/1000;

%average radius : 2.2886/2/1000 m , 2.9582*10^3 drops per m^3
probability=0.081; %0.081

drop=[0 0 0 0 0 0];
%drop : [r v A m(unit : kg) F h (coalenscnece num)]
%        1 2 3 4            5 6        7

% now=c+drop_count;
life=0; %added
ancestor=0;

times_gone=0;
times=0;
while 1
    times_gone=times*dt;
    times=times+1;
    ask_create=rand();
    if drop_count==0 || ask_create<=probability
        k=rand();
        func=@(r) p(r)-k;
        d_rand=fzero(func,1.5)/1000; % unit : m
        if d_rand<0
            d_rand=fzero(func,3)/1000;
        end

        drop_count=drop_count+1;
        
        drop(1,:,drop_count)=[d_rand/2 0 pi*(d_rand/2)^2 4/3*pi*(d_rand/2)^3*1000 4/3*pi*(d_rand/2)^3*1000*g h0];
        life(drop_count)=1; %added
    end
    
    for i = 1:drop_count

        r_temp=drop(life(i),1,i);
        v_temp=drop(life(i),2,i);
        a_temp=drop(life(i),3,i);
        m_temp=drop(life(i),4,i);
        F_temp=drop(life(i),5,i);
        h_temp=drop(life(i),6,i);
        
        v_temp=(F_temp*dt+drop(life(i),4,i)*v_temp)/m_temp;
        h_temp=h_temp-v_temp*dt;
        rho_a=-0.0977*h_temp/1000+1.225;
        F_temp=m_temp*g-1/2*Cd*rho_a*a_temp*v_temp^2;%-rho_a*4/3*pi*r_temp^3*g;%-rho_a*4/3*a_temp*drop(nth,1,1)*g; % add buyoancy
        
        life(i)=life(i)+1; %added
        drop(life(i),:,i)=[r_temp v_temp a_temp m_temp F_temp h_temp];
        
        
        if h_temp<0 %h<0 : arrived on land, remove all the data of the drop
            ancestor=ancestor+1
%             if drop(1,1,i)>=1.10/1000 && drop(1,1,i)<=1.20/1000
                range=(times-life(i)+1)*dt:dt:(times)*dt;
                subplot(2,1,1)
                plot(range,drop(1:life(i),5,i))
                hold on
                subplot(2,1,2)
                plot(range,drop(1:life(i),1,i))
                hold on
%                 range=0:dt:(life(i)-1)*dt;
%                 plot(range,drop(1:life(i),2,i))
%                 hold on
%             end
            drop_count=drop_count-1;
            
            for j = i:drop_count
                drop(:,:,j)=drop(:,:,j+1);
                life(j)=life(j+1); %added
                for k = life(j)+1:length(drop(:,1,1))
                    drop(k,:,j)=[0 0 0 0 0 0];
                end
            end
        end
        
    end
    
    
    for i = 1:drop_count-1 % check if they can do coalescence ?
        h1=drop(life(i),6,i); h2=drop(life(i+1),6,i+1); r1=drop(life(i),1,i); r2=drop(life(i+1),1,i+1); %a
        
        if (h2-h1)<abs(r2+r1) || h2-h1<0 % if 2 drops are in adequate coalescence range
            temp=drop(life(i+1),:,i+1); % save last phase of drops that will be removed
            determine_h=min(h1,h2);
            %after coalescence; -> determine mass, velocity... etc
            r_final=(temp(1)^3+drop(life(i),1,i)^3)^(1/3);
            m_final=temp(4)+drop(life(i),4,i);
            v_final=(temp(4)*temp(2)+drop(life(i),4,i)*drop(life(i),2,i))/m_final;
            rho=-0.0977*determine_h/1000+1.225;
            
                
            
            if r_final>=r_max  %split! because it's too big to be stable
                k=rand();
                func=@(r) p(r)-k;
                d_rand=fzero(func,1.5)/1000; % unit : m
                if d_rand<=0
                    d_rand=fzero(func,3)/1000;
                end
%                 r_temp=drop(1,1,i);
                r_temp=d_rand/2;
                F_temp=4/3*pi*r_temp^3*1000*g-1/2*Cd*rho*pi*r_temp^2*v_final^2;%-rho*4/3*r_temp^3*pi*g;
                
                drop(life(i),:,i)=[r_temp v_final pi*(r_temp)^2 4/3*pi*(r_temp)^3*1000 F_temp determine_h];
                
            else %normal case/ used conservation of momentum...
                a_final=r_final^2*pi;
                F_final=m_final*g-1/2*Cd*rho*a_final*v_final^2;%-rho*4/3*r_final^3*pi*g; %add buyoancy
                drop(life(i),:,i)=[r_final v_final pi*r_final^2 m_final F_final determine_h];
            end
            
            for j= i+1:drop_count-1
                drop(:,:,j)=drop(:,:,j+1);
                life(j)=life(j+1);
            end
            drop_count=drop_count-1; 
        end
    end
    
    for l=1:drop_count %destroy weird data or config
        if abs(drop(life(l),2,l))>15
            drop(life(l),5,l)=0;
%             for j=l:drop_count-1
%                 drop(:,:,j)=drop(:,:,j+1);
%                 life(j)=life(j+1);
%             end
%             drop_count=drop_count-1;
        end
    end
        
end