n=10;
% city=rand(n,2); % random init of coordinates of cities
                % 'city' matrix and 'n' can be provided explicitly instead 
                % of random init of coordinates of cities
city =[0.0640    0.4192;
    0.5628    0.9510;
    0.8773    0.9927;
    0.3052    0.1991;
    0.6523    0.8703;
    0.2235    0.6982;
    0.0350    0.2255;
    0.6443    0.6001;
    0.5395    0.0192;
    0.9961    0.7132];
                
% city =[0.5851    0.9890;
% 0.5749    0.8064;
% 0.2065    0.0430;
% 0.0571    0.5424;
% 0.2295    0.2647];

%order = [1 2 3 5 4];
%order = [3 7 5 1 2 6 8 4];
order = [7 1 6 2 5 3 10 8 9 4];

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
plot(city(9,1),city(9,2),'r*','MarkerSize',12);
plot(city(10,1),city(10,2),'g*','MarkerSize',12);
xlabel('x-coordinate')
ylabel('y-coordinate')
legend('1','2','3','4','5', '6', '7', '8', '9', '10');
%%

% plot(complex(cos(a*0),sin(a*0)),'ko','MarkerSize',4);
% plot(complex(cos(a*1),sin(a*1)),'ko','MarkerSize',4);
% plot(complex(cos(a*2),sin(a*2)),'ko','MarkerSize',4);
% plot(complex(cos(a*3),sin(a*3)),'ko','MarkerSize',4);
% plot(complex(cos(a*4),sin(a*4)),'ko','MarkerSize',4);
% plot(complex(cos(a*5),sin(a*5)),'ko','MarkerSize',4);
% plot(complex(cos(a*6),sin(a*6)),'ko','MarkerSize',4);
% plot(complex(cos(a*7),sin(a*7)),'ko','MarkerSize',4);
% plot(complex(cos(a*8),sin(a*8)),'ko','MarkerSize',4);
% plot(complex(cos(a*9),sin(a*9)),'ko','MarkerSize',4);