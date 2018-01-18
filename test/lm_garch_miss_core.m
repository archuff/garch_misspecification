function stat = lm_garch_miss_core(T,h,param,vart,varx,sres,nblag)

% Purpose: 
%   Auxiliary function which has to be used with lm_garch_miss.m

dof = size(vart,2);
v = zeros(T,dof);
x = zeros(T,nblag);
dlnB = zeros(T,dof);
dvv = zeros(dof);
dxv = zeros(nblag,dof);
dxx = zeros(nblag,nblag);

for t = 2:T,
    temp = zeros(dof,t);
    temp2 = zeros(nblag,t);
    for j = 1:t-1,
        temp(:,j) = param(3)^(j-1)*vart(t-j+1,:)';
        temp2(:,j) = param(3)^(j-1)*varx(t-j,:)';
    end
    v(t,1:dof) = (1/h(t)).*sum(temp,2)';
    x(t,:) = (1/h(t)).*sum(temp2,2)';
    dlnB(t,:) = (sres(t)*(v(t,:)));
    dvv = dvv + (v(t,:)'*v(t,:));
    dxv = dxv + (x(t,:)'*v(t,:));
    dxx = dxx + (x(t,:)'*x(t,:));
end

sdlnB = sum(dlnB,1);
V = dvv-(dxv'*inv(dxx)*dxv);
stat = (1/2)*(sdlnB*inv(V)*sdlnB');

