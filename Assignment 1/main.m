
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
T = 500; % in ms
del_t = 0.1;
I_c = 2.7e-9; % calculated in part 1 (a)
m = (T/del_t); % 500/0.1 + 1;
n = 10;
I = ones(n,1);
for i = 1:n
    I(i,1) = (1 + i*alpha)*I_c;
end
A = I.*ones(1,m);
V_t = runge_kutta(A, V_thresh, El, del_t);
t = 0:0.1:500;
%disp(size(V_t));
plot(V_t(2,:), 'Linewidth', 2);
 