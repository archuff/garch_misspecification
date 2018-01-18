function stat = ho_test_core(T,h,param,vart,varx,sres)

v  = zeros(T,2);
dlnB = zeros(2,T);
x  = zeros(T,3);
dvv = zeros(2,2);
dxv = zeros(3,2);
dxx = zeros(3,3);

for t = 2:T,
    temp = zeros(2,t);
    temp2 = zeros(3,t);   
    for j = 1:t-1,
        temp(:,j) = param(3)^(j-1)*vart(t-j,:)';
        temp2(:,j) = param(3)^(j-1)*varx(t-j,:)';
    end
    v(t,:) = (1/h(t)).*sum(temp,2);
    dlnB(:,t) = sres(t)*v(t,:);
    x(t,:) = (1/h(t)).*sum(temp2,2)';
    dvv = dvv + (v(t,:)'*v(t,:));
    dxv = dxv + (x(t,:)'*v(t,:));
    dxx = dxx + (x(t,:)'*x(t,:));
end

sdlnB = sum(dlnB,2)';
V = (1/T)*((sres'*sres)/T)*(dvv-dxv'*inv(dxx)*dxv);
stat = (1/T)*(sdlnB*(V)^(-1)*sdlnB');

