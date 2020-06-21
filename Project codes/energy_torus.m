function [L,grad_L]=energy_torus(cities,R, r, theta_z,k,l)
  n=size(cities,1);
 
  A=5; B = 5; C = 20; D = 2; F = 5;   % parameters
  
  a=0;b=0;c=0;d=0; f=0;
  z0=complex(0,1); 
  theta = angle(theta_z);
  coord_x = (R + r.*cos(theta(:,1))).*cos(theta(:,2));
  coord_y = (R + r.*cos(theta(:,1))).*sin(theta(:,2));
  coord_z = r.*sin(theta(:,1));
  for i=1:n
    a=a + ((R - sqrt(coord_x(i).^2 + coord_y(i).^2)).^2 + coord_z(i).^2 - r^2);
    z1=theta_z(i,2); %theta_z(i,2)
    b=b+abs((z1^n)-1)^2;
    %f = f + ((abs(theta_z(i,1))^2)-1)^2 + ((abs(theta_z(i,2))^2)-1)^2;
    for j=1:n
      z2=theta_z(j,2);
      e=z2*(z1');
      d=d+dist(cities(i,:),cities(j,:))*exp((-imag(sqrt(e))^2)/k);
      c=c+(abs(z1-z2))^2;
    end
  end
  if l==1
    disp('Value of terms:');
    disp(a);
    disp(b);
    disp(c);
    disp(d);
    disp(f);
  end
  L=A*a+B*b -C*c + D*d + F*f;    % tour energy
 
  grad_L=zeros(n,2);            % initializing gradient vector
      
  %theta_1
  for i=1:n
      grad_L(i,2)= 2*F*((abs(theta_z(i,2))^2)-1)*theta_z(i,2) + B*((n*z0*imag((theta_z(i,2)/abs(theta_z(i,2)))^n))/(theta_z(i,2)')) - 2*A*complex(sin(theta(i,2)), cos(theta(i,2)))*((R - sqrt(coord_x(i).^2 + coord_y(i).^2))*(coord_x(i)*sin(theta(i,2)) - coord_y(i)*cos(theta(i,2)))/(sqrt(coord_x(i)^2 + coord_y(i)^2)));
      grad_L(i,1) = 2*F*((abs(theta_z(i,1))^2)-1)*theta_z(i,1) - 2*A*complex(sin(theta(i,1)), cos(theta(i,1)))*((R - sqrt(coord_x(i).^2 + coord_y(i).^2))*(coord_x(i)*sin(theta(i,1)*cos(theta(i,2)) + coord_y(i)*sin(theta(i,1))*sin(theta(i,2)) )/(sqrt(coord_x(i)^2 + coord_y(i)^2))) + coord_z(i)*cos(theta(i,1)));
      for j=1:n
         e=(theta_z(j,2)*(theta_z(i,2)'))/(abs(theta_z(j,2))*abs(theta_z(i,2)));
         grad_L(i, 2)= grad_L(i, 2)+(D*z0*dist(cities(i,:),cities(j,:))*exp((-imag(sqrt(e))^2)/k)*imag(e))/(4*k*(theta_z(i,2)'));
         q=imag(e');
         grad_L(i, 2)= grad_L(i, 2)-((C*z0*q)/(theta_z(i,2)'));
      end

  end
end
