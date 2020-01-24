function[V_t] = runge_kutta(A, n, m, del_t)
b = -30e-9/300e-12;
c = 1/300e-12;
El = -70e-3;
d = b*El;
h = del_t;
V_t = zeros(n, m+1);
V_thresh = 20e-3;
%Assuming zero initial condition
%solving V'(t) = b*V(t) + c*I(t) + d
for i = 1:m
 k1 = b.*V_t(:,i) + c.*A(:,i) + d;
 y1 = V_t(:,i) + k1*h/2;
 k2 = b.*y1 + c.*A(:,i) + d; %Assuming current same for period of h/2
 V_t(:,i+1) = y1 + k2*h/2;  
 if V_t(:,i+1) >= V_thresh 
    V_t(:,i+1) = El;
 end
end
V_t = V_t(:,(2:m+1));
end