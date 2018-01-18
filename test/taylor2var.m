function [taylorResults] = taylor2var(x1,x2,k)

% Auxilary function returning a matrix where vectors are the Taylor expansion 
% components of a function with more than one argument
% Author: Thomas Chuffart
% Mail: thomas.chuffart@univ-fcomte.fr

    [lx1,cx1] = size(x1);
    [lx2,cx2] = size(x2);

    tot = 0;
    for i = 1:k,
        m = cx1 + i - 1; 
        tot  = combmj(m,i) + tot ;
    end
    
    taylorResults = zeros(lx1,tot);
    count = 1; 
    for i = 1:cx1,
        itaylor = 1;
        if itaylor == k+1, break, end         
        temp1 = x1(:,i); 
        taylorResults(:,count) = temp1;
        count = count +1;
        for j = i:cx2,
            itaylor = 2;
            if itaylor == k+1, break, end 
            temp2 = temp1.*x2(:,j);
            taylorResults(:,count) = temp2;
            count = count+1;
            for l = j:cx2,
                itaylor = 3;
                if itaylor == k+1, break, end                 
                temp3 = temp2.*x2(:,l);
                taylorResults(:,count) = temp3;
                count = count+1;
                for m = l:cx2,
                    itaylor = 4;
                    if itaylor == k+1, break, end 
                    temp4 = temp3.*x2(:,m);
                    taylorResults(:,count) = temp4;
                    count = count+1; 
                    for n = m:cx2,
                        itaylor = 5;
                        if itaylor == k+1, break, end               
                        temp5 = temp4.*x2(:,n);
                        taylorResults(:,count) = temp5;
                        count = count+1; 
                        for o = n:cx2,
                            itaylor = 6;
                            if itaylor == k+1, break, end               
                            temp6 = temp5.*x2(:,o);
                            taylorResults(:,count) = temp6;
                            count = count+1;                          
                        end
                    end
                end
            end       
        end
    end       
end

