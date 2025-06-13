 %% Main --- Non-stoch targeted exploration

run initialize_nonstoch.m
% noise energy bound
% ws=normrnd(mu_w, sqrt(1), nx, T);
% gammaw=0;
% for k=1:T
%     gammaw=gammaw+(ws(:,k)'*ws(:,k));
% end

D0tilde=1e1*eye(nphi);
D0=kron(D0tilde,eye(nx));
D0inv=inv(D0);
    
DTtilde=1e0*eye(nphi);
DT=kron(DTtilde,eye(nx));
DThalf=chol(DT);
sims=[];
%%
utilde=200*ones(nu,L);
u=utilde;
gammaw=5;
run newprior_nonstoch.m
% gammaw=1000;
run lmi_test.m

%%
% gammaw=ceil(gammaw);
run exploration_nonstoch.m
% run exploration_nonstoch_YALMIP.m
% disp(u);
un=u/norm(u);
nsn=[0,gammaw, un];
sims=[sims;nsn];
disp("done");
%%
Uc=zeros(L,T);
tset=0:T-1;
for i=1:L
    cosi=cos(2*pi*freqs(i)*tset);
    Uc(i,:)=u(i)*cosi;
end
Un=sum(Uc);

wstoch=[];
Xn=zeros(nx,T+1);
for n=1:T
    ws=[-sqrt(gammaw/T)*cos(Xn(1,n));0;0;0];
    wstoch=[wstoch,ws];
%     xn=A*Xn(:,n)+B*Un(:,n)+ws(:,n);
    xn=A*Xn(:,n)+B*Un(:,n)+ws;
%     xn=A*Xn(:,n)+B*Un(:,n);
    Xn(:,n+1)=xn;
end

dn=zeros(nphi);
for k=1:T
    dn=dn+[Xn(:,k);Un(:,k)]*[Xn(:,k);Un(:,k)]';
end

Pninv=kron(dn,eye(nx));

tmp=zeros(20,1);
for k=1:T
    tmp=tmp+kron([Xn(:,k);Un(:,k)],eye(nx))*Xn(:,k+1);
end
theta_n=Pninv\tmp;
ndist=(thetatr-theta_n)'*(thetatr-theta_n);
disp(max(eig(ndist)));
% pn=reshape(theta_n,[nx,nphi]);
% An=pn(:,1:nx);
% Bn=pn(:,nx+1);

xns=0;
for k=1:T
    xns=xns+(Xn(:,k+1)'*Xn(:,k+1));
end


G=gammaw+(theta_n'*(Pninv)*theta_n)-xns;
% save("ga.mat","sims","-append");

%%
sig=10;
ut=200*ones(nu,L);
%%
run exploration_stoch.m

% disp(us);
usn=us/norm(us);
stoch=[1,sig, usn];
sims=[sims;stoch];
% save("ga.mat","sims","-append");

%%
cdel=chi2inv(1-0.01,(nx^2)+(nx*nu));
us=(us/norm(us))*norm(u);
Uc=zeros(L,T);
for i=1:L
    cosi=cos(2*pi*freqs(i)*tset);
    Uc(i,:)=us(i)*cosi;
end
Us=sum(Uc);

Xs=zeros(nx,T+1);

for n=1:T
%     xs=A*Xs(:,n)+B*Us(:,n)+ws(:,n);
    xs=A*Xs(:,n)+B*Us(:,n)+wstoch(:,n);
    % xs=A*Xs(:,n)+B*Us(:,n);
    Xs(:,n+1)=xs;
end

ds=zeros(nphi);
for k=1:T
    ds=ds+[Xs(:,k);Us(:,k)]*[Xs(:,k);Us(:,k)]';
end

Psinv=kron(ds,eye(nx));

tmp=zeros(20,1);
for k=1:T
    tmp=tmp+kron([Xs(:,k);Us(:,k)],eye(nx))*Xs(:,k+1);
end

theta_s=Psinv\tmp;
sdist=((thetatr-theta_s)'*(thetatr-theta_s));
disp(max(eig(sdist)));
%%
gammaw=100;
sig=gammaw/T;
ut=1000*ones(nu,L);

%%
gtemp=inf;
for ct=1:1
    run exploration_nonstoch.m
    if (abs(gtemp-gammae)/gammae)<1e-2
        break;
    end
    gtemp=gammae;
end

