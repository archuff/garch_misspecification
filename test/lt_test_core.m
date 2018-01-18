function stat = lt_test_core(T,h,param,vart,varx,sres,nblag)

dof = 2;
v = vart;
dvv = zeros(dof);
dxv = zeros(nblag,dof);
dxx = zeros(nblag,nblag);
x = zeros(T,nblag);
dlnB = zeros(T,dof);

for t = 2:T,
        temp = zeros(nblag,t);
        for j = 1:t-1,
            temp(:,j) = param(3)^(j-1)*varx(t-j,:)';
        end
        x(t,:) = 1/h(t).*sum(temp,2);
        v(t,:) = 1/h(t).*vart(t-1,:);
        dlnB(t,:) = (sres(t)*(v(t,:)));
        dvv = dvv + (v(t,:)'*v(t,:));
        dxv = dxv + (x(t,:)'*v(t,:));
        dxx = dxx + (x(t,:)'*x(t,:));
end  

sdlnB = sum(dlnB,1)';
V = dvv-dxv'*inv(dxx)*dxv;
k = sum(sres(2:T,1).^2)/(T-1);
V = k*V;
stat = sdlnB'*inv(V)*sdlnB;
