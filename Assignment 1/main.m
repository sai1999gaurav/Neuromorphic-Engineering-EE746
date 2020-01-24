
%% %Problem 1 (b)
A=input('Input Matrix :\n');
[n,m] = size(A);
T = input('total time of consideration : \n');
V_t = runge_kutta(A, n, m, T/m);
disp(V_t);
%% Problem 1(c)
alpha = 0.1;
T = 500; % in ms
del_t = 0.1;
I_c = 2.7e-9; % calculated in part 1 (a)
m = T/del_t; % 500/0.1;
n = 4;
I = ones(n,1);
for i = 1:n
    I(i,1) = (1 + 2*i*alpha)*I_c;
end
A = zeros(n, m+1);
A(:,(2:m+1)) = I.*ones(1,m);
V_t = runge_kutta(A, n, m, del_t);
t = 0:0.1:500;
 for i = 1:n
 plot(t, V_t(i,:), 'Linewidth', 2);
 hold on;
 end
 
