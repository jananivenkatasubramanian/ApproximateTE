%% Get prior estimates A0 and B0 from data
% count=1;
% while count>0
%     theta0=(mvnrnd(thetatr,D0inv))';
%     tr=(thetatr-theta0)'*D0*(thetatr-theta0);
%     if tr<10
%         count=0;
%     end
% end

theta0=thetatr+0.5e-2;
% tr=(thetatr-theta0)'*D0*(thetatr-theta0);

p0=reshape(theta0,[nx,nphi]);
A0=p0(:,1:nx);
B0=p0(:,nx+1);
% disp(p0);


%% transfer matrices to compare with the non-stochastic exploration
Vh=[];

for i=1:nphi
    vt=inv(exp(1i*2*pi*freqs(i))*eye(nx)-A0)*B0;
    vi=[vt;eye(nu)];
    Vh=[Vh,vi];
end
