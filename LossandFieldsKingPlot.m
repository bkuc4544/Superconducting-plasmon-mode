clear all;
clc;


%Constants
%King and Wiltse: Surface-Wave Propagation


T = 4.2;                        % in K
Tc = 9.26;                      % in K
n = 5.56e28;                    % in /m3, Ref: Mermin, Ashcroft, Solid State Physics, pg 5
nn = n*(T/Tc)^4;
ns = n - nn;

freq = [2:105]*10^9; %Frequency 
lambda_L = 32e-9;               % Penetration Depth
c = 3e8;
a = 0.08382*1e-3;
mu0 = 4*pi*10^-7;
tau = 1e-12;
e = 1.6e-19;
m = 9.1e-31;
ns = m/(mu0 * lambda_L^2 * e^2);
sigma1 = nn*(e^2)*tau./(m*(1+(2*pi*freq*tau).^2));                                                          % For SC Nb
sigma2 = ns * e^2 ./ (m * 2*pi*freq) + nn*(e^2)*tau^2*(2*pi*freq).^2./(m*(2*pi*freq+tau^2*(2^3*pi^3*freq.^3)));    % For SC Nb
%sigmac = 7e6;                                           % For normal Nb
%sigmac = 5.8e7;                                         % For normal Cu
%sigmac = 3.8e7;                                         % For normal Al
%sigmac = 1.67e7;                                        % For normal Zn
muc = mu0;
%kc=sqrt(2*pi*freq*muc*sigmac)*exp(-1i*pi/4);            % For normal Nb/Cu
sigmac = sigma1 - 1i*sigma2;                             % For SC Nb
kc=sqrt(-1i*2*pi*muc.*freq.*sigmac);                         % For SC Nb
gammac = kc;

for pp=1:length(freq)
    w1 = 2*pi*freq(pp);
    k = w1/c;

    f = @(g) abs(besselh(0,((g(1)+1i*g(2))*a))./besselh(1,(g(1)+1i*g(2))*a)-1i*k.^2./((g(1)+1i*g(2))*kc(pp)));  % The parameterized function.
    X(pp,:) = fminsearch(@(g) f(g),[15, 15]);

    gr(pp)=X(pp,1);
    gi(pp)=X(pp,2);
end
    
%% Plot Loss

gamma=gr+1i*gi;
k=2*pi*freq/c;
h=sqrt(k.^2-gamma.^2);

theta=angle(gamma*a);

beta=abs(gamma.*a).^2./(2*k*a.^2).*sin(2*theta); %Loss in Np/m

figure(1)
plot(freq/10^9,-beta*8.685889638, 'Linewidth',2)             
xlabel('Freq (GHz)')
ylabel('Loss (dB/m)')
title('Loss vs Frequency for Nb')
set(gca,'Fontsize',18)
grid on

%% Plot fields
A=1;
Nf=49; % Choose Frequency
r = linspace(a, 10000*a, 100000);
Hphi=1i*A*k(Nf)^2/(2*pi*freq(Nf)*mu0*gamma(Nf))*besselh(1,(gr(Nf)+1i*gi(Nf))*r);
Er=1i*A*h(Nf)/(gamma(Nf))*besselh(1,(gr(Nf)+1i*gi(Nf))*r);
Ez=A*besselh(0,(gr(Nf)+1i*gi(Nf))*r);


figure(2)
subplot(2,2,1),plot(r*100,abs(Hphi),r*100,abs(Er)/10^3,'Linewidth',2)
title(['Frequency = ' num2str(freq(Nf)/10^9) ' GHz'])
xlabel('Distance (cm)')
ylabel('Amplitude (kV/m or A/m)')
set(gca,'Fontsize',18)
legend('H_\phi','E_r')
grid on

subplot(2,2,2),plot(r*100,abs(Ez),'Linewidth',2)
title(['Frequency = ' num2str(freq(Nf)/10^9) ' GHz'])
xlabel('Distance (cm)')
ylabel('Amplitude (V/m)')
set(gca,'Fontsize',18)
legend('E_z')
grid on

subplot(2,2,3),plot(r*100,abs(Hphi)/max(abs(Hphi)),r*100,abs(Er)/max(abs(Er)),'--',r*100,abs(Ez)/max(abs(Ez)),'Linewidth',2)
%title(['Frequency = ' num2str(ff(Nf)/10^9) ' GHz'])
xlabel('Distance (cm)')
ylabel('Norm. Amplitude')
set(gca,'Fontsize',18)
legend('H_\phi','E_r','Ez')
grid on

subplot(2,2,4)
loglog(r*100,abs(Hphi),r*100,abs(Er)/10^3,r*100,abs(Ez)/10^3, 'Linewidth',2)
%title(['Frequency = ' num2str(ff(Nf)/10^9) ' GHz'])
xlabel('Distance (cm)')
ylabel('Amplitude (kV/m or A/m)')
set(gca,'Fontsize',18)
legend('H_\phi','E_r','E_z')
grid on


figure(3)
Zs = Ez./Hphi;
plot(r*100, real(Zs), r*100, imag(Zs))
title(['Frequency = ' num2str(freq(Nf)/10^9) ' GHz'])
xlabel('Radial Distance (cm)')
ylabel('Surface Impedance (\Omega)')
set(gca,'Fontsize',18)
grid on

Frequency=freq(Nf);

%save Fields Frequency r Hphi Er Ez % in GHz meters A/m V/m V/m

%%
load BareWireLareRadius
figure(3)
rhfss=BareWireLargeRadius(:,1)*100;
loglog(rhfss,BareWireLargeRadius(:,2:end),rhfss,BareWireLargeRadiusMode2(:,2:end))

Escale=90;

loglog(r*100,abs(Ez)/Escale,rhfss,BareWireLargeRadius(:,[2 4])/10^10,rhfss,BareWireLargeRadiusMode2(:,3)/10^10)
loglog(r*100,abs(Ez)/Escale,'k',rhfss,BareWireLargeRadiusMode2(:,3)/10^10,'r--','Linewidth',2)
xlim([1 400])


figure(4)
loglog(r*100,abs(Er)/Escale,rhfss,BareWireLargeRadius(:,[8 10])/10^10,rhfss,BareWireLargeRadiusMode2(:,9)/10^10)
loglog(r*100,abs(Er)/Escale,'k',rhfss,BareWireLargeRadiusMode2(:,9)/10^10,'r--','Linewidth',2)
xlim([1 400])

figure(5)
loglog(r,abs(Ez)/Escale,'k',rhfss/100,BareWireLargeRadiusMode2(:,3)/10^10,'r--',...
    r,abs(Er)/Escale,'b',rhfss/100,BareWireLargeRadiusMode2(:,9)/10^10,'g--','Linewidth',2)
xlim([0.01 4])

legend('E_z Theory','E_z FEM','E_r Theory','E_r FEM') 
xlabel('Distance (m)')
ylabel('E-field (V/m)')
set(gca,'Fontsize',18)


