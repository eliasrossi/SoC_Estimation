function [xh,zh] = func_EKF(eqdif,funch,x0,z,u,par)

% eqdif - state space equation
% funch - output equation
% x0 - initial conditions
% z - measured data
% u - input data
% par - parameters
% mat0 - initialization of covariance matrices

n=length(z);
dt=par.dt;
nx=par.nx;
nz=par.nz;
nu=par.nu;

A=zeros(nx,nx);
B=zeros(nx,nu);
C=zeros(nz,nx);

dh=1e-9;

P0 = par.P0;
Q = par.Q;
R = par.R;
P=P0;

xh = zeros(nx,n); % the first collum is the initial condition
zh = zeros(nz,n);
xh(:,1) = x0; % initial condition
% 
for j=1:n-1
    if abs(det(P))>det(P0) || det(P)==0 || isnan(det(P)) 
        P=P0;
    end

    for i=1:nx
        xh_delta=xh(:,j);
        xh_delta(i)=xh_delta(i)+dh;
        A(:,i)=(eqdif(xh_delta,u(:,j),par)-eqdif(xh(:,j),u(:,j),par))/dh;
    end
    
    for i=1:1:nu
        u_delta=u(:,j);
        u_delta(i)=u_delta(i)+dh;
        B(:,i)=(eqdif(xh(:,j),u_delta,par)-eqdif(xh(:,j),u(:,j),par))/dh;
    end
    
    for i=1:nx
        xh_delta=xh(:,j);
        xh_delta(i)=xh_delta(i)+dh;
        C(:,i)=(funch(xh_delta,u(:,j),par)-funch(xh(:,j),u(:,j),par))/dh;
    end
    F=A;
    H=C;
    %% STATE COVARIANCE COMPUTATION
    % state prediction covariance
    Pb = F*P*F' + Q; 
    if abs(det(Pb))>det(P0) || det(Pb)==0 || isnan(det(Pb)) 
        Pb=P0;
    end
    % innovation covariance
    S = R + H*Pb*H';
%     if abs(det(S))>1e10 || det(S)==0 || isnan(det(S)) 
%         S=R;
%     end
    if det(S)==0 || isnan(det(S)) 
        S=R;
    end
%     det(S)

    % Filter gain
    W = Pb*H'/S;
    %% ESTIMATION OF THE STATE
    % state prediction
    xb = eqdif(xh(:,j),u(:,j),par);
    % measurement prediction
    zb=funch(xb,u(:,j),par);
    % measurement residual 
    ni = z(:,j+1) - zb;
    % updated state estimate
    xplus = xb + W*ni;
    
    xplus(xplus < par.LB) = par.LB(xplus < par.LB); % bounds
    xplus(xplus > par.UB) = par.UB(xplus > par.UB); % bounds

    xh(:,j+1)=xplus;
    zh(:,j+1)=funch(xh(:,j+1),u(:,j+1),par);
        %% Updated state covariance    
    P = Pb - W*S*W';
end

end