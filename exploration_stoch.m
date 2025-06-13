%% Solve targeted exploration

cvx_begin SDP
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

variable us(nu,L)

variable gammae(1)

% estimates A0 and B0
% X_T and Phi_T to be generated in loop here
% L1 and L2 to be generated in loop with random multiplier first

Ue=blkdiag(us(:,1),us(:,2),us(:,3),us(:,4),us(:,5));
Util=blkdiag(ut(:,1),ut(:,2),ut(:,3),ut(:,4),ut(:,5));
%% Exploration energy

Se_stoch_energy=[gammae, ones(1,nphi)*(Ue'); Ue*ones(nphi,1), gammae*eye(nu*nphi)];
% Se_energy=[gammae, u1; u1, gammae];

%% Exploration LMI

% exploration block matrices ----------
Sexp1_stoch=Vh*(Ue*Util'+Util*Ue'-Util*Util')*Vh';
Sexp_stoch=Sexp1_stoch-(DTtilde*L*sig*chi2inv(1-0.01,(nx^2)+(nx*nu))/(T));
%%
tol=1e-1;

minimize gammae
% subject to
Se_stoch_energy>=0;
% Se_energy>=-tol*eye(nphi+1);
% gammae>=0;
% gammae>=-tol;
Sexp_stoch>=0;

cvx_end

% disp('cond(S_exp)='); disp(cond(S_exp));
% disp('min(eig(S_exp))='); disp(min(eig(S_exp)));
% disp('max(eig(S_exp))='); disp(max(eig(S_exp)));


ut=us;

% Ut=Ue;


%%





