function [V_t] = runge_kutta(I, V_thresh, El, del_t)
[n,m] = size(I);
%solving CV'(t) = -gl(V(t)-El) + I
gl = 30e-9;
C = 300e-12;
h = del_t;
V_t = zeros(n, m);   %Initial V = El
V_t(:,1) = V_t(:,1) + El.*ones(n,1);
k1=ones(n,1);
k2=ones(n,1);
for i = 1:m-1
    %Runge-Kutta-Algo
k1=(-gl*(V_t(:,i)-El*ones(n,1))+I(:,i))/C;
k2=(-gl*((V_t(:,i)+k1.*h)-El*ones(n,1))+I(:,i+1))/C;
V_t(:,i+1) = V_t(:,i) + ((k1+k2)/2).*h;
%V_t_i = h*(I(:,i-1) - gl*(V_t(:,i-1) - El))/C + V_t(:,i-1);
%V_t(:,i) = h*(((I(:,i-1) - gl*(V_t(:,i-1) - El))/C + (I(:,i) - gl*(V_t_i - El))/C)/2) + V_t(:,i-1);
V_t(:,i+1) = (V_t(:,i+1) >= V_thresh)*El + (V_t(:,i+1) < V_thresh).*V_t(:,i+1); % Compare with Threshold 
end
end