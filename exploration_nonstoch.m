%% Solve targeted exploration

cvx_begin SDP quiet
cvx_solver SDPT3
% cvx_precision high

%% VARIABLES
%Dt lower bound
% variable Dt(nphi,nphi)
% Dt == semidefinite(nphi);

% variable u1(nu,1) 
% variable u2(nu,1) 
% variable u3(nu,1) 
% variable u4(nu,1) 
% variable u5(nu,1)

variable u(nu,L)

variable gammae(1)

% estimates A0 and B0
% X_T and Phi_T to be generated in loop here
% L1 and L2 to be generated in loop with random multiplier first

Ue=blkdiag(u(:,1),u(:,2),u(:,3),u(:,4),u(:,5));

%% Exploration energy

Se_energy=[gammae, ones(1,nphi)*(Ue'); Ue*ones(nphi,1), gammae*eye(nu*nphi)];
% Se_energy=[gammae, u1; u1, gammae];

%% Exploration LMI

% exploration block matrices ----------
Xhat=zeros(1,T*nx); %\hat{X}^\top
Phihat=zeros(nphi,T); %\hat{\Phi}
Yhat=zeros((nx*nphi)+(nx*nphi)^2,T*nx*nx*nphi); %\hat{Y}
L21=zeros(1,T*nx);
L1=zeros(nphi,T); %L_1 for convex relaxation of \Phi
L2=zeros((nx*nphi)+(nx*nphi)^2,T*nx*nx*nphi); % L_2 for convex relaxation of Y
xi=[];xi1=[];
phi=[];phi1=[];
xtemp=x0;xtemp1=x0;

for i=1:L
    for j=1:T
        xit=A0*xtemp+B0*u(:,i)*Uc(i,j);
        xit1=A0*xtemp1+B0*utilde(:,i)*Uc(i,j);
        xi=[xi,xit'];
        xi1=[xi1,xit1'];
        phi=[phi,[xit;u(:,i)*Uc(i,j)]];
        phi1=[phi1,[xit1;utilde(:,i)*Uc(i,j)]];
        xtemp=xit;xtemp1=xit1;
    end
    xtemp=x0;xtemp1=x0;
    Xhat=Xhat+xi; %this is \hat{X}^\top summed over frequencies
    Phihat=Phihat+phi; %this is \hat{\Phi} summed over frequencies
    L1=L1+phi1; %L_1 for convex relaxation of \Phi
    L21=L21+xi1; % used to define L_2
    Yhat=Yhat+[DThalf'*kron(Xhat,eye(nx*nphi));kron(kron(Phihat,eye(nx)),eye(nx*nphi))];
    L2=L2+[DThalf'*kron(L21,eye(nx*nphi));kron(kron(L1,eye(nx)),eye(nx*nphi))];
    phi=[];phi1=[];
    xi=[];xi1=[];
end

Se1_1=kron((Phihat*L1'+L1*Phihat'-L1*L1'),eye(nx))-gammaw*DT;
Se1_2=zeros(nx*nx*nphi*nphi);
Se1=blkdiag(Se1_1,Se1_2);

Se2=Yhat*L2'+L2*Yhat'-L2*L2';

S_exp=Se1+Se2;

%%
tol=1e-1;

minimize gammae
% subject to
Se_energy>=0;
% Se_energy>=-tol*eye(nphi+1);
% gammae>=0;
% gammae>=-tol;
S_exp>=0;

cvx_end

% disp('cond(S_exp)='); disp(cond(S_exp));
% disp('min(eig(S_exp))='); disp(min(eig(S_exp)));
% disp('max(eig(S_exp))='); disp(max(eig(S_exp)));


utilde=u;

% Ut=Ue;


%%





