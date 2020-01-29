%% %Problem 1 (b)
V_thresh = 20e-3;
El = -70e-3;
I=input('Input Current Matrix :\n');
T = input('total time of consideration : \n');
V_t = runge_kutta(I, V_thresh, El, T/size(I,2));
disp(V_t);
%% Problem 1(c)
V_thresh = 20e-3;
El = -70e-3;
alpha = 0.1;
T = 500e-3;
del_t = 0.1e-3;
I_c = 2.7e-9; % calculated in part 1 (a)
m = (T/del_t); % 500/0.1 + 1;
n = 10;
I = ones(n,1);
%I(:,1)=zeros(n,1);
for i = 2:n
    I(i,1) = (1 + i*alpha).*I_c;
end
A = I*ones(1,m);
V_t = runge_kutta(A, V_thresh, El, del_t);
t = 0:0.1:T;
%disp(size(V_t));
hold on
plot(V_t(2,:), 'Linewidth', 2)
plot(V_t(4,:), 'Linewidth', 2)
plot(V_t(6,:), 'Linewidth', 2)
plot(V_t(8,:), 'Linewidth', 2)
legend('n=2','n=4','n=6','n=8')
title('LIF Model')
hold off
%% Problem 2(c)
El = -70e-3;
alpha = 0.1;
del_t = 0.1e-3;
I_c = 2.7e-9; % calculated in part 1 (a)
kz = 0;
C = 0;
Er = 0;
Et = 0;
a = 0;
b = 0;
c= 0;
d = 0;
U_steady = 0;
V_steady = 0;
V_thresh = 0;
n = 10;
I_app = 400e-12;
T = 100e-3;
m = (T/del_t); % 500/0.1 + 1;
type1 = 1;
type2 = 2;
type3 = 3;
I = I_app*ones(n,1);
A = I*ones(1,m);
V_t1 = runge_kutta_fourth_order(A, type1, del_t);
V_t2 = runge_kutta_fourth_order(A, type2, del_t);
V_t3 = runge_kutta_fourth_order(A, type3, del_t);
t = 0:0.1:T;
%disp(size(V_t));
hold on
plot(V_t1(4,:), 'Linewidth', 2)
plot(V_t2(4,:), 'Linewidth', 2)
plot(V_t3(4,:), 'Linewidth', 2)
%plot(V_t(4,:), 'Linewidth', 2)
%plot(V_t(6,:), 'Linewidth', 2)
%plot(V_t(8,:), 'Linewidth', 2)
legend('RS','IB','CH')
title('Izhikevich Model')
hold off