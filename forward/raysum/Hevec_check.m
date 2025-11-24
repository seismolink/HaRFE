function [phase1,mult_ml,errflag]=evec_check(evec,invec,phase1)
% Consistency check for eigenvectors. Checks that evec (current
% eigenvectors) is consistent with invec (eigenvector for current
% phase calculated at previous step). phase1 (indicating which
% column in evec corresponds to invec) may be modified if a 90-degree
% rotation has occurred; if the sign is wrong, mult is set to -1.

mult_ml=1;
errflag=false;

%       index: 1 if downgoing, 4 if up
index_ml=fix(((abs(1.).*sign(real(phase1)-3.5+eps))./2. + 0.5).*3.+1.);
%          echk=evec(1:3,index:(index+2))  [relevant eigenvectors]
echk = evec(1:3,index_ml:(index_ml+2));
%          comp_list=real(invec'*echk) [check which evec is closest]
comp_list = real(invec'*echk);
[maxval_ml,maxindex] = max(abs(comp_list));
%       check if the phase matches
if(maxindex ~=(phase1-index_ml+1))
    phase1=maxindex+index_ml-1;
end
%       check if the sign matches
if (real(comp_list(maxindex)) < 0)
    mult_ml=-1;
end
%       check if the evecs don't line up. This shouldn't happen.
if(maxval_ml < 0.99)
%     disp('ERROR')
    errflag=true;
end

end