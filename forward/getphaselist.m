function [ray] = getphaselist(Nlay,noreflect)

% Computation of the raypathes and their wavetype
% -----------------------------------------------
all_perm = 1;
ray_path_base = (Nlay+1)-(1:1:Nlay);
ray(1).path = ray_path_base;
ray(1).type = ones(size(ray_path_base));
indx=2;
for j = 1:Nlay
    for i = 1:Nlay
        if j == 1 && i == 1
            Nrp = length(ray_path_base);
            A = [2 3]'*ones(1,Nrp-1);
            B=nan(2^size(A,2),size(A,2));
            for k=1:size(A,2)
                B(:,k)=repmat(reshape(repmat(A(:,k)',2^(size(A,2)-k),1),[],1),2^(k-1),1);
            end
            Nrp_perm = size(B);
            Nperm = Nrp_perm(1);
            for l=1:Nperm
                for k = 1:Nlay-1
                    i1 = k+1;
                    i2 = Nlay;
                    ray(indx).path=ray_path_base;
%                     ray(indx).type=[1 B(l,:)];
                    ray(indx).type = ones(size(ray_path_base));
                    ray(indx).type(i1:i2) = B(l,(i1:i2)-1);
                    indx=indx+1;
                end
            end
        end
        if j >1 && (j+i) < (Nlay+2)
            if noreflect(Nlay+2-j-i) || noreflect(Nlay-i+1)
                continue
            end
            ref_vec_1 = ((Nlay+2-j-i):(Nlay-i));
            ref_vec_2 = abs(-(Nlay-i):-(Nlay+2-j-i));
            ref_vec=[ref_vec_1 ref_vec_2];
            raypath=ray_path_base(1:(i+j-1));
            raypath((i+j):(i+j+length(ref_vec)-1))=ref_vec;
            raypath((i+j+length(ref_vec)):(length(ref_vec)+length(ray_path_base)))=ray_path_base((i+j):end);
            if all_perm
            Nrp = (length(ref_vec)+length(ray_path_base));
            A = [1 2 3]'*ones(1,Nrp);
            B=nan(3^size(A,2),size(A,2));
            for k=1:size(A,2);
                B(:,k)=repmat(reshape(repmat(A(:,k)',3^(size(A,2)-k),1),[],1),3^(k-1),1);
            end
            B(:,(i+j):length(ref_vec_1))=B(:,(i+j):length(ref_vec_1))+3;
            else
                Nrp = (length(ref_vec)+length(ray_path_base))-i-j+1;
                A = [2 3]'*ones(1,Nrp);
                B=nan(2^size(A,2),size(A,2));
                for k=1:size(A,2);
                    B(:,k)=repmat(reshape(repmat(A(:,k)',2^(size(A,2)-k),1),[],1),2^(k-1),1);
                end
                
                Nrp_perm = size(B);
                B(:,(i+j):(Nrp_perm(2)+(i+j)-1))=B;
                B(:,1:(i+j-1))=ones(Nrp_perm(1),(i+j-1));
            end
                B(:,(i+j):(length(ref_vec_1)+i+j-1))=B(:,(i+j):(length(ref_vec_1)+i+j-1))+3;
       %     end
            Nrp_perm = size(B);
            Nperm = Nrp_perm(1);
            for l=1:Nperm
                if B(l,1) == 1 && B(l,Nrp_perm(2)) ~= 1
            ray(indx).path=raypath;
            ray(indx).type=B(l,:);
            indx=indx+1;
                end
            end
        end
    end
end
% zu den strukturvariablen path und type kommen im Laufe der Berechnung
% noch time und amp_coe hinzu

end