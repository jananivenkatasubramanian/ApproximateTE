%% Initialization - Sequential learning and control
clear all;
clc;

% system order
nx=4;
nu=1;
nphi=nx+nu;

x0=zeros(nx,1);

% A, B and frequencies
freqs=[0,0.1,0.2,0.3,0.4];
L=size(freqs,2);
A=[0.49 0.49 0 0;0 0.49 0.49 0; 0 0 0.49 0.49;0 0 0 0.49];
B=[0; 0; 0; 0.49];
thetatr=[A(:);B];

% u and w parameters to obtain prior
mu_u=0;
sigma_u=sqrt(100);

mu_w=0;
sigma_w=1;

% number of trials, number of rollouts
N_r=100;

% Exploration time
T=100;

% Omega_T set of possible frequencies
Omega_T=[];
for i=1:T
    Omega_T(i)=(i-1)/T;
end

% u_k^{(i)}
Uc=zeros(L,T);
tset=0:T-1;
for i=1:L
    cosi=cos(2*pi*freqs(i)*tset);
    Uc(i,:)=cosi;  
end


w1=zeros(nx,T-2);
% w2=(sqrt(10)/2)*ones(nx,1);
w2=[sqrt(10);0;0;0];
ws=[w1,w2,[0;0;0;0]];
wn=ws;
gammaw=0;
for k=1:T
    gammaw=gammaw+(ws(:,k)'*ws(:,k));
end

%%
% gammatest=0;
% for k=1:T
%     gammatest=gammatest+(wstoch(:,k)'*wstoch(:,k));
% end