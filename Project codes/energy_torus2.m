function [L,grad_L]=energy_torus2(cities,R, r, theta_z,k,l)
  n=size(cities,1);
 
  A= 0; B = 4; C = 1e-2; D = 7; F = 0.5;   % parameters A = 0.1, F = 0.5
  
  a=0;b=0;c=0;d=0; f=0;
  z0=complex(0,1); 
  theta = angle(theta_z);
  coord_x = (R + r.*cos(theta(:,1))).*cos(theta(:,2));
  coord_y = (R + r.*cos(theta(:,1))).*sin(theta(:,2));
  coord_z = r.*sin(theta(:,1));
  sp = zeros(n,n);
  for i=1:n
    a=a + ((R - sqrt(coord_x(i).^2 + coord_y(i).^2)).^2 + coord_z(i).^2 - r^2);
    z1=theta_z(i,2)/abs(theta_z(i,2)); %theta_z(i,2)
    b=b+abs((z1^n)-1)^2;
    f = f + ((abs(theta_z(i,1))^2)-1)^2 + ((abs(theta_z(i,2))^2)-1)^2;
    for j=1:n
      z2=theta_z(j,2)/abs(theta_z(j,2));
      e=z2*(z1');
      d=d+dist(cities(i,:),cities(j,:))*exp((-imag(sqrt(e))^2)/k);
      %c=c+(abs(z1-z2))^2;
      if (i ~= j)
          sp(i,j) = 2*pi*sqrt( r^2*(theta(i,1) - theta(j,1))^2 + (theta(i,2) - theta(j,2))^2*(R + r*cos(theta(i,1)))*(R + r*cos(theta(j,1))));
          c = c + sp(i,j); 
      end
    end
  end
  if l==1
%     disp('Value of terms:');
%     disp(a);
%     disp(b);
%     disp(c);
%     disp(d);
%     disp(f);
  end
  L=A*a+B*b -(C/2)*c + (D/2)*d + F*f;    % tour energy
 
  grad_L=zeros(n,2);            % initializing gradient vector
      
  %theta_1
  for i=1:n
      %grad_L(i,2)= 2*F*((abs(theta_z(i,2))^2)-1)*theta_z(i,2) + B*((n*z0*imag((theta_z(i,2)/abs(theta_z(i,2)))^n))/(theta_z(i,2)')) - 2*A*complex(sin(theta(i,2)), cos(theta(i,2)))*((R - sqrt(coord_x(i).^2 + coord_y(i).^2))*(coord_x(i)*sin(theta(i,2)) - coord_y(i)*cos(theta(i,2)))/(sqrt(coord_x(i)^2 + coord_y(i)^2)));
      grad_L(i,2) = 2*F*((abs(theta_z(i,2))^2)-1)*theta_z(i,2) + B*((n*z0*imag((theta_z(i,2)/abs(theta_z(i,2)))^n))/(theta_z(i,2)')) - 2*A*complex(sin(theta(i,2)), cos(theta(i,2)))*((R - sqrt(coord_x(i).^2 + coord_y(i).^2))*(coord_x(i)*sin(theta(i,2)) - coord_y(i)*cos(theta(i,2)))/(sqrt(coord_x(i)^2 + coord_y(i)^2)));     
      %grad_L(i,1) = 2*F*((abs(theta_z(i,1))^2)-1)*theta_z(i,1) - 2*A*complex(sin(theta(i,1)), cos(theta(i,1)))*((R - sqrt(coord_x(i).^2 + coord_y(i).^2))*(coord_x(i)*sin(theta(i,1)*cos(theta(i,2)) + coord_y(i)*sin(theta(i,1))*sin(theta(i,2)) )/(sqrt(coord_x(i)^2 + coord_y(i)^2))) + coord_z(i)*cos(theta(i,1)));
      grad_L(i,1) = 2*F*((abs(theta_z(i,1))^2)-1)*theta_z(i,1) - 2*A*complex(sin(theta(i,1)), cos(theta(i,1)))*((R - sqrt(coord_x(i).^2 + coord_y(i).^2))*(coord_x(i)*sin(theta(i,1)*cos(theta(i,2)) + coord_y(i)*sin(theta(i,1))*sin(theta(i,2)) )/(sqrt(coord_x(i)^2 + coord_y(i)^2))) + coord_z(i)*cos(theta(i,1)));
      for j=1:n
         e=(theta_z(j,2)*(theta_z(i,2)'))/(abs(theta_z(j,2))*abs(theta_z(i,2)));
         grad_L(i, 2)= grad_L(i, 2)+(D*z0*dist(cities(i,:),cities(j,:))*exp((-imag(sqrt(e))^2)/k)*imag(e))/(4*k*(theta_z(i,2)'));
         %q=imag(e');
         if (i ~= j)
             grad_L(i, 2)= grad_L(i, 2) - C*((2*pi)^2/sp(i,j))*((R + r*cos(theta(i,1)))*(R + r*cos(theta(j,1)))*(theta(i,2)*complex(sin(theta(i,2)), cos(theta(i,2))) - theta(j,2)*complex(sin(theta(j,2)), cos(theta(j,2))) ) );
         end
      end

  end
end
