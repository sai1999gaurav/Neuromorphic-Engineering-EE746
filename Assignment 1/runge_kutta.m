function[V_t] = runge_kutta(I, V_thresh, El, del_t)
[n,m] = size(I);
%solving CV'(t) = -gl(V(t)-El) + I
gl = 30e-9;
C = 300e-12;
h = del_t;
V_t = zeros(n, m) + El.*ones(n,1);  %Initial V = El
for i = 2:m
    %Runge-Kutta-Algo
    V_t_i = h*(I(:,i-1) - gl*(V_t(:,i-1) - El))/C + V_t(:,i-1);
    V_t(:,i) = h*(((I(:,i-1) - gl*(V_t(:,i-1) - El))/C + (I(:,i) - gl*(V_t_i - El))/C)/2) + V_t(:,i-1);
    for j = 1:n
     if (V_t(j,i)>=V_thresh)
      V_t(j,i) = El;
     end
    end
end 

end