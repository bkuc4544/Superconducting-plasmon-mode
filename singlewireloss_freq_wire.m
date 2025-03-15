clear all;
%close all;
clc;

%THz Sommerfeld wave propagation on a single metal wire

%H = besselh(nu,K,Z) computes the Hankel function H^(K)_?(z)

%Constants;

sigma_c=5.8*10^7; %copper
a=5e-6;
fv=[5:5:100]*10^9;
mu0 = 4*pi*10^-7;
c = 3e-8;
for pp=1:length(fv)
f=fv(pp);
omega=2*pi*f;
kc=sqrt(omega*mu0*sigma_c)*exp(-1i*pi/4);
wl=c./f;
k=2*pi./wl;
g=[15:0.1:30]*exp(1i*pi*0.3627);

LH=besselh(0,1,g*a)./besselh(1,1,g*a);
RH=1i*k.^2./(g*kc);

%%
plot(abs(g),real(RH),'b',abs(g),real(LH),'r',...
    abs(g),imag(RH),'k--',abs(g),imag(LH),'r--');

% 9.0725 +19.7124i

f = @(g) abs(besselh(0,1,((g(1)+1i*g(2))*a))./besselh(1,1,(g(1)+1i*g(2))*a)-1i*k.^2./((g(1)+1i*g(2))*kc));  % The parameterized function.
                        % The parameter.
X(pp,:) = fminsearch(@(g) f(g),[15 15]);

h(pp)=sqrt(-(X(pp,1)+1i*X(pp,2)).^2+k^2);

end



plot(fv/10^9,-imag(h),'Linewidth',2)
xlabel('Frequency (GHz)')
ylabel('Amplitude Absorption (m^{-1})')


plot(fv/10^9,20*log10(exp(imag(h))),'Linewidth',2)
xlabel('Frequency (GHz)')
ylabel('Attenuation (dB/m)')
title(['Single Wire Radius ' num2str(a*100) 'cm'])
return
%%
N=10;
r=a:a/10:20*a;
plot(r*100,abs(besselh(0,1,((X(N,1)+1i*X(N,2))*r)))/abs(besselh(0,1,((X(N,1)+1i*X(N,2))*a))),'Linewidth',2)
xlabel('Distance (cm)')
ylabel('Normalized Transverse Field')
title([num2str(fv(N)/10^9) ' GHz Single Wire Radius ' num2str(a*100) 'cm'])
xlim([0 max(r)*100])
