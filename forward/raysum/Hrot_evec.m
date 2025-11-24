function [evec]=rot_evec(evec,r)
%       Rotate degenerate post-rotation eigenvectors into a consistent
%       coordinate system

x3(1)=complex(0.,0.);
x3(2)=complex(0.,0.);
x3(3)=complex(1.,0.);

for i=1:6
%             prod_i = R*evec(1:3,i)
    prod_ml(:,i) = r*evec(1:3,i);
end

%       Check if SH is off-horizontal in rotated frame.
%       If (R*u3).x3 > 0 or (R*u6).x3 > 0...
if ((abs(dot(prod_ml(:,3),x3)) > eps) || (abs(dot(prod_ml(:,6),x3)) > eps))

%         Downgoing set:
theta=atan2(real(dot(prod_ml(:,3),x3)), real(dot(prod_ml(:,2),x3)));

a(1,1)=complex(cos(theta));
a(1,2)=complex(-sin(theta));
a(2,1)=complex(sin(theta));
a(2,2)=complex(cos(theta));
%          evec(:,2:3)*A
evec2 = evec(:,2:3)*a;
evec(:,2:3) = evec2;

%         Upgoing set:
theta=atan2(real(dot(prod_ml(:,6),x3)), real(dot(prod_ml(:,5),x3)));
a(1,1)=complex(cos(theta));
a(1,2)=complex(-sin(theta));
a(2,1)=complex(sin(theta));
a(2,2)=complex(cos(theta));
%          evec(:,5:6)*A
evec2 = evec(:,5:6)*a;
evec(:,5:6) = evec2;

end

end