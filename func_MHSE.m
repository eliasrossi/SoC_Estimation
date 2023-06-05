function [xh,zh,xest] = func_MHSE(eqdif,funch,x0,z,u,par)

% eqdif - state equation (linear or nonlinear) of the system
% funch - output equation (linear or nonlinear) of the system

nx=par.nx;
n=par.n;
nz=par.nz;
N=par.N; % tamanho da janela
xh=zeros(nx,n);
zh = zeros(nz,n);
xest=zeros(nx,n-N);

xh(:,1)=x0;

for j=2:N
    xh(:,j)=eqdif(xh(:,j-1),u(:,j-1),par);
    zh(:,j)=funch(xh(:,j),u(:,j),par); 
end

xp=x0;

for j=N+1:1:n
    j
%     mincost = @(x) costfcn(eqdif,funch,x,xp,z(:,j-N:j),par,u(:,j-N:j));
    mincost = @(x) costfcn2(eqdif,funch,x,xp,z(:,j-N:j),par,u(:,j-N:j));
    options=optimoptions('lsqnonlin','OptimalityTolerance',1e-8,'UseParallel',true,'Display','Off');
    xh_tN=lsqnonlin(mincost,xp,par.LB,par.UB,options);
    xest(:,j-N)=xh_tN;
    xp=eqdif(xh_tN,u(:,j-N),par);
    xh_aux=xh_tN;
    for i=N:-1:1
        xh_aux=eqdif(xh_aux,u(:,j-i),par);
    end
    xh(:,j)=xh_aux;
    zh(:,j)=funch(xh(:,j),u(:,j),par);
end

end