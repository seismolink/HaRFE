function [eval,evec]=reorder_evec(evalr,evali,evecr,eveci,index_ml)

for j=1:6
    eval(j)=complex(evalr(index_ml(j)),evali(index_ml(j)));
    for k=1:6
        evec(k,j)=complex(evecr(k,index_ml(j)),eveci(k,index_ml(j)));
    end
end

end