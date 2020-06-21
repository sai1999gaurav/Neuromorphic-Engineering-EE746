function [L,grad_L]=energy_3(cities,z,k,C,l)
  n=size(cities,1);
 
  A=0.5;B=7;D=10;   % parameters
  
  a=0;b=0;c=0;d=0;
  z0=complex(0,1); 
  
  for i=1:n,
    a=a+((abs(z(i))^2)-1)^2;
    z1=z(i)/abs(z(i));
    c=c+abs((z1^n)-1)^2;
    for j=1:n,
      z2=z(j)/abs(z(j));
      e=z2*(z1');
      b=b+dist(cities(i,:),cities(j,:))*exp((-imag(sqrt(e))^2)/k);
      d=d+(abs(z1-z2))^2;
    end
  end
  if l==1
    disp('Value of terms:');
    disp(a);
    disp(b);
    disp(c);
    disp(d);
  end
  L=A*a+(B/2)*b+C*c-(D/2)*d;    % tour energy
 
  grad_L=zeros(n,1);            % initializing gradient vector
      
  for i=1:n
      grad_L(i)=2*A*((abs(z(i))^2)-1)*z(i)+C*((n*z0*imag((z(i)/abs(z(i)))^n))/(z(i)'));
      for j=1:n
         e=(z(j)*(z(i)'))/(abs(z(j))*abs(z(i)));
         grad_L(i)=grad_L(i)+(B*z0*dist(cities(i,:),cities(j,:))*exp((-imag(sqrt(e))^2)/k)*imag(e))/(4*k*(z(i)'));
         q=imag(e');
         grad_L(i)=grad_L(i)-((D*z0*q)/(z(i)'));
      end

  end
end
