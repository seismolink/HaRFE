function [cr]=rot_tensor(cc,r)
% Rotate tensor CC according to rotator matrix R. Result is output
% as tensor CR.

cr = zeros(3,3,3,3);

for i=1:3
for j=1:3
for k=1:3
for l=1:3

for a=1:3
for b=1:3
for c=1:3
for d=1:3

cr(i,j,k,l)=cr(i,j,k,l) + r(a,i).*r(b,j).*r(c,k).*r(d,l).* cc(a,b,c,d);

end
end
end
end

end
end
end
end

end