%% Initialize coordinates of cities
%n=input('number of cities:');
clc;
close all;
clear all;


n=8;
% city=rand(n,2); % random init of coordinates of cities
                % 'city' matrix and 'n' can be provided explicitly instead 
                % of random init of coordinates of cities
% city =[0.0640    0.4192;
%     0.5628    0.9510;
%     0.8773    0.9927;
%     0.3052    0.1991;
%     0.6523    0.8703;
%     0.2235    0.6982;
%     0.0350    0.2255;
%     0.6443    0.6001;
%     0.5395    0.0192;
%     0.9961    0.7132];
                
% city =[0.5851    0.9890;
% 0.5749    0.8064;
% 0.2065    0.0430;
% 0.0571    0.5424;
% 0.2295    0.2647];

%order = [1 2 3 5 4];
order = [3 7 5 1 2 6 8 4];
%order = [7 1 6 2 5 3 10 8 9 4];

% city =[0.2458    0.5678;
%     0.6762    0.6130;
%     0.2613    0.8070;
%     0.9615    0.8808;
%     0.4976    0.1600;
%     0.4278    0.0892];

city =[0.7641    0.6762;
     0.8    0.6058;
     0.1253    0.863;
     0.3645    0.2920;
     0.5235    0.8300;
     0.9708    0.4900;
     0.4387    0.7674;
     0.8119    0.0133];
       % 0.8961    0.3485];                

l=0;
a=(2*pi)/n;
figure(1);
hold on;
plot(city(1,1),city(1,2),'bx','MarkerSize',12);
plot(city(2,1),city(2,2),'kx','MarkerSize',12);
plot(city(3,1),city(3,2),'rx','MarkerSize',12);
plot(city(4,1),city(4,2),'gx','MarkerSize',12);
plot(city(5,1),city(5,2),'mx','MarkerSize',12);
plot(city(6,1),city(6,2),'yx','MarkerSize',12);
plot(city(7,1),city(7,2),'b*','MarkerSize',12);
plot(city(8,1),city(8,2),'k*','MarkerSize',12);
% plot(city(9,1),city(9,2),'r*','MarkerSize',12);
% plot(city(10,1),city(10,2),'g*','MarkerSize',12);
xlabel('x-coordinate')
ylabel('y-coordinate')
legend('1','2','3','4','5', '6', '7', '8' ); %, '9', '10');
%title('Coordinates of cities');
%hold off;
corr_sum = 0;
fault_lim = 1;
corr_tot = 100;
%%
for iter = 1:corr_tot

runs=1;
L_hist=zeros(runs,1);

ap=1;
Cv=0.4;     % Value of parameter C
kv=0.1;     % Value of parameter k

res=zeros(size(kv,2),ap);
for w=1:ap,
    sp=zeros(size(kv,2),1);
    min_E=zeros(size(kv,2),1);
    for e=1:size(kv,2),
        for k=1:runs,
        %disp(k);
        % Parameters initialization
        mag=2*rand(n,1);
        theta=2*pi*rand(n,1);
        z=complex(mag.*cos(theta),mag.*sin(theta));  % initializing z 
        
        z_init=z;
        iter=12000;   % total no. of iterations  12000
        fac=0.1;     % Noise Decay factor
        alpha=0.01;  % delta t
        ann = 1 - 1e-2;
        L_energy=zeros(iter,1);
        z_e=zeros(iter,n);
        gr_L = ones(iter, 1);


        %% Update of zj
        iter1 = 0;
        i = 0;
        while (max(abs(gr_L)) > 1e-2)
            i = i + 1;
            iter1 = iter1 + 1;
            tn=2*pi*rand(n,1);
            zn=complex(cos(fac*tn),sin(fac*tn));
            fac=ann*fac;
            z=z.*zn;
            [L,gr_L]=energy_3(city,z,kv(e),Cv,l);
            L_energy(i)=L;
            for j=1:n
                z(j)=z(j)-alpha*gr_L(j);
                z_e(i,j)=z(j);
            end
        end
        
        disp('iterations');
        disp(iter1);
        
        %% Plotting evolution of energy L
        if (k<=2)
%             figure(2);
%             hold on;
%             plot(L_energy,'color',rand(1,3));
%             title('Evolution of L');
%             xlabel('Iterations');
%             ylabel('Energy L');
        end
        %disp('Final L');
        %disp(L);

        if (mod(k,1)==0)
            %disp('gr');
            %disp(gr_L);
        end
        L_hist(k)=round(L*10000)/(10000);
        
        end

    E=unique(L_hist);
    count=zeros(size(E,1),1);
    t=0;
    for i=1:size(E,1)
        if isnan(E(i))==0
            for j=1:size(L_hist,1),
                if L_hist(j)==E(i)
                    count(i)=count(i)+1;
                    t=t+1;
                end
            end
        end
    end
    %disp([E,count])
    %disp(t)     % total no.of non NaN final energy
    ind=find(isnan(E)==0);
    %figure(4);
    %bar(E(ind),count(ind)/t);
    %disp([E(ind),count(ind)]);
    %xlabel('Energy');
    %ylabel('Probability');

    i=find(E==min(E));
    %disp('E=');
    %disp(E);
    %disp('min E=');
    %disp(min(E));
    min_E(e)=min(E);
    %sp(e)=count(i)/t;
    end

% figure(5);
% plot(kv,sp,'color',rand(1,3));
% disp('The success prob is');
% disp('C=');disp(Cv);
% disp('kv  sp min-E=');
% disp([kv',sp,min_E]);
% hold on;
% axis([0.4,2,0,1]);
% xlabel('Value of Parameter');ylabel('Success Probability');
res(:,w)=sp;
end

%% Plotting initial and final z on unit circle 
l=1;
[L,gr_L]=energy_3(city,z,kv(e),Cv,l);

% figure(4);
% plot(real(z_e(1,1)),imag(z_e(1,1)), '*');
% hold on
% plot(real(z_e(1,:)),imag(z_e(1,:)));
% disp('total iterations');
% disp(iter1);
        
figure(3);
x=0:0.01*pi:2*pi;
z1=complex(cos(x),sin(x));
plot(z1,'c:');
hold on;
axis([-1.5,1.5,-1.5,1.5]);

%% To view evolution of nodes 
% pause(5);
% for k=1:500
%     %disp(k);
%     plot(z1,'c:');
%     plot(z_e(k,1),'bx','MarkerSize',12);
%     plot(z_e(k,2),'kx','MarkerSize',12);
%     plot(z_e(k,3),'rx','MarkerSize',12);
%     plot(z_e(k,4),'gx','MarkerSize',12);
%     plot(z_e(k,5),'mx','MarkerSize',12);
% 
%     %plot(z_e(k,6),'bo','MarkerSize',6);
%     %plot(z_e(k,7),'ko','MarkerSize',6);
%     %plot(z_e(k,8),'ro','MarkerSize',6);
%     %plot(z_e(k,9),'go','MarkerSize',6);
%     %plot(z_e(k,10),'mo','MarkerSize',6);
% 
%     %plot(z_e(k,11),'b*','MarkerSize',10);
%     %plot(z_e(k,12),'k*','MarkerSize',10);
%     %plot(z_e(k,13),'r*','MarkerSize',10);
%     %plot(z_e(k,14),'g*','MarkerSize',10);
%     %plot(z_e(k,15),'m*','MarkerSize',10);
% 
%     
% %     plot(z_init(1),'bo','MarkerSize',6);
% %     plot(z_init(2),'ko','MarkerSize',6);
% %     plot(z_init(3),'ro','MarkerSize',6);
% %     plot(z_init(4),'go','MarkerSize',6);
% %     plot(z_init(5),'mo','MarkerSize',6);
% %     
%     pause(1e-3);
%     plot(z_e(k,1),'wx','MarkerSize',12);
%     plot(z_e(k,2),'wx','MarkerSize',12);
%     plot(z_e(k,3),'wx','MarkerSize',12);
%     plot(z_e(k,4),'wx','MarkerSize',12);
%     plot(z_e(k,5),'wx','MarkerSize',12);
% 
%     %plot(z_e(k,6),'wo','MarkerSize',6);
%     %plot(z_e(k,7),'wo','MarkerSize',6);
%     %plot(z_e(k,8),'wo','MarkerSize',6);
%     %plot(z_e(k,9),'wo','MarkerSize',6);
%     %plot(z_e(k,10),'wo','MarkerSize',6);
% 
%     %plot(z_e(k,11),'w*','MarkerSize',10);
%     %plot(z_e(k,12),'w*','MarkerSize',10);
%     %plot(z_e(k,13),'w*','MarkerSize',10);
%     %plot(z_e(k,14),'w*','MarkerSize',10);
%     %plot(z_e(k,15),'w*','MarkerSize',10);
% 
% end




% plot(real(z_e(end,6)),imag(z_e(end,6)),'yx','MarkerSize',14);
% 
% plot(real(z_e(end,7)),imag(z_e(end,7)),'b*','MarkerSize',14);
% plot(real(z_e(end,8)),imag(z_e(end,8)),'k*','MarkerSize',14);
% plot(real(z_e(end,9)),imag(z_e(end,9)),'r*','MarkerSize',14);
% plot(real(z_e(end,10)),imag(z_e(end,10)),'g*','MarkerSize',14);

%legend('0', '1','2','3','4','5', '6', '7', '8');
%disp(angle(theta_z(:,2)));
[sort_ang, sort_ind] = sort(angle(z));
%disp(sort_ang);
%disp(sort_ind);
b_ind = find(order == sort_ind(1));
b_ind = rem(b_ind, n) + 1;
corr = 1;
fault = 0;
for i = 2:n
    if (sort_ind(i) ~= order(b_ind))
        fault = fault + 1;
        %break        
    end
    if (fault > fault_lim)
        corr = 0;
       break 
    end
    b_ind = rem(b_ind, n) + 1;   
end

disp('corr');
disp(corr);
corr_sum = corr_sum + corr;
if (corr == 1)
    plot(real(z(1)),imag(z(1)),'bx','MarkerSize',14);
plot(real(z(2)),imag(z(2)),'kx','MarkerSize',14);
plot(real(z(3)),imag(z(3)),'rx','MarkerSize',14);
plot(real(z(4)),imag(z(4)),'gx','MarkerSize',14);
plot(real(z(5)),imag(z(5)),'mx','MarkerSize',14);
plot(real(z(6)),imag(z(6)),'yx','MarkerSize',14);
plot(real(z(7)),imag(z(7)),'r*','MarkerSize',14);
plot(real(z(8)),imag(z(8)),'g*','MarkerSize',14);
    break
end
end

disp('success prob');
disp(corr_sum/corr_tot);

%% Roots of unity
plot(complex(cos(a*0),sin(a*0)),'ko','MarkerSize',4);
plot(complex(cos(a*1),sin(a*1)),'ko','MarkerSize',4);
plot(complex(cos(a*2),sin(a*2)),'ko','MarkerSize',4);
plot(complex(cos(a*3),sin(a*3)),'ko','MarkerSize',4);
plot(complex(cos(a*4),sin(a*4)),'ko','MarkerSize',4);
plot(complex(cos(a*5),sin(a*5)),'ko','MarkerSize',4);
plot(complex(cos(a*6),sin(a*6)),'ko','MarkerSize',4);
plot(complex(cos(a*7),sin(a*7)),'ko','MarkerSize',4);
% plot(complex(cos(a*8),sin(a*8)),'ko','MarkerSize',4);
% plot(complex(cos(a*9),sin(a*9)),'ko','MarkerSize',4);
