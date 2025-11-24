function [eval,evec]=sort_evec(evalr,evali,evecr,eveci)

% Got burned by using a DATA statement here...
nrp=0;
nrn=0;
nip=0;
nin=0;
realpos = zeros(1,6);
realneg = zeros(1,6);
imagpos = zeros(1,6);
imagneg = zeros(1,6);
index_ml = zeros(1,6);

% Divide eigenvalues up into real positive, real negative, complex
% positive, complex negative.
for i=1:6
    %         Is it real?
    if(abs(evali(i)./evalr(i)) < eps*100)
        evali(i)=0;
        for j=1:6
            eveci(j,i)=0;
        end
        if (evalr(i) >= 0.)
            nrp=nrp+1;
            realpos(nrp)=i;
        else
            nrn=nrn+1;
            realneg(nrn)=i;
        end
    else
        if (evali(i) >= 0.)
            nip=nip+1;
            imagpos(nip)=i;
        else
            nin=nin+1;
            imagneg(nin)=i;
        end
    end
end

realpos(realpos==0) = [];
realneg(realneg==0) = [];
imagpos(imagpos==0) = [];
imagneg(imagneg==0) = [];
%        Sort sub-groups
% [realpos]=mysort(evalr,realpos,nrp);
% [realneg]=mysort(evalr,realneg,nrn);
% [imagpos]=mysort(evali,imagpos,nip);
% [imagneg]=mysort(evali,imagneg,nin);
if ~nrp
else
[~,I]=sort(evalr(realpos));
realpos = realpos(I);
end
if ~nrn
else
[~,I]=sort(evalr(realneg));
realneg = realneg(I);
end
if ~nip
else
[~,I]=sort(evali(imagpos));
imagpos = imagpos(I);
end
if ~nin
else
[~,I]=sort(evali(imagneg));
imagneg = imagneg(I);
end

%        Assemble sub-ranges
for i=1:nip
    index_ml(i)=imagpos(i);
end
for i=1:nrp
    index_ml(i+nip)=realpos(i);
end
for i=1:nin
    index_ml(i+nip+nrp)=imagneg(nin-i+1);
end
for i=1:nrn
    index_ml(i+nip+nrp+nin)=realneg(nrn-i+1);
end


[eval,evec]=Hreorder_evec(evalr,evali,evecr,eveci,index_ml);

end