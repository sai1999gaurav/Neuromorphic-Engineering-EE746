function [V_t] = runge_kutta_fourth_order(I, type, del_t)
[n,m] = size(I);
kz = 0.7e-6*(type-2)*(type-3)/2-1.2e-6*(type-1)*(type-3)+1.5e-6*(type-1)*(type-2)/2;
C = 100e-12*(type-2)*(type-3)/2-150e-12*(type-1)*(type-3)+50e-12*(type-1)*(type-2)/2;
Er = (-60e-3*(type-2)*(type-3)/2+75e-3*(type-1)*(type-3)-60e-3*(type-1)*(type-2)/2)*ones(n,1);
Et= (-40e-3*(type-2)*(type-3)/2+45e-3*(type-1)*(type-3)-40e-3*(type-1)*(type-2)/2)*ones(n,1);
a = 30*(type-2)*(type-3)/2-10*(type-1)*(type-3)+30*(type-1)*(type-2)/2;
b = (-2*(type-2)*(type-3)/2-5*(type-1)*(type-3)+1*(type-1)*(type-2))*1e-9/2;
c = ((-50*(type-2)*(type-3)/2+56*(type-1)*(type-3)-40*(type-1)*(type-2))*1e-3/2);
d = ((100*(type-2)*(type-3)/2-130*(type-1)*(type-3)+150*(type-1)*(type-2))*1e-12/2);
V_thresh = (35*(type-2)*(type-3)/2-50*(type-1)*(type-3)+25*(type-1)*(type-2))*1e-3/2;
U_steady = (-34.4*(type-2)*(type-3)/2-170.85*(type-1)*(type-3)+15.66*(type-1)*(type-2))*1e-12/2;
V_steady = (-42.8*(type-2)*(type-3)/2+40.83*(type-1)*(type-3)-44.33*(type-1)*(type-2))*1e-3/2;

h = del_t;
V_t = zeros(n, m);   %Initial V = El
U_t = zeros(n,m);
U_t(:,1) = U_t(:,1) + U_steady.*ones(n,1);
V_t(:,1) = V_t(:,1) + V_steady.*ones(n,1);
k1=ones(n,1);
k2=ones(n,1);
k3=ones(n,1);
k4=ones(n,1);
k5=ones(n,1);
k6=ones(n,1);
k7=ones(n,1);
k8=ones(n,1);
for i = 1:m-1
%Runge-Kutta-Algo
k1=((kz*(V_t(:,i)-Er).*(V_t(:,i)-Et)-U_t(:,i)+I(:,i))/C)*h;
k2=((kz*(V_t(:,i)+k1/2-Er).*(V_t(:,i)+k1/2-Et)-U_t(:,i)+I(:,i))/C)*h;
k3=((kz*(V_t(:,i)+k2/2-Er).*(V_t(:,i)+k2/2-Et)-U_t(:,i)+I(:,i))/C)*h;
k4=((kz*(V_t(:,i)+k3-Er).*(V_t(:,i)+k3-Et)-U_t(:,i)+I(:,i))/C)*h;
k5=(a*(b*(V_t(:,i)-Er)-U_t(:,i)))*h;
k6=(a*(b*(V_t(:,i)-Er)-U_t(:,i)-k5/2))*h;
k7=(a*(b*(V_t(:,i)-Er)-U_t(:,i)-k6/2))*h;
k8=(a*(b*(V_t(:,i)-Er)-U_t(:,i)-k7))*h;
V_t(:,i+1) = V_t(:,i) + ((k1+2*k2+2*k3+k4)/6);
U_t(:,i+1) = U_t(:,i) + ((k5+2*k6+2*k7+k8)/6);

%V_t_i = h*(I(:,i-1) - gl*(V_t(:,i-1) - El))/C + V_t(:,i-1);
    %V_t(:,i) = h*(((I(:,i-1) - gl*(V_t(:,i-1) - El))/C + (I(:,i) - gl*(V_t_i - El))/C)/2) + V_t(:,i-1);
    for j = 1:n
       % || (V_t(j,i)==El 
     if (V_t(j,i+1)>=V_thresh) 
         V_t(j,i+1) = c;
         U_t(j,i+1) = U_t(j,i+1) + d;
     end
end 
end