%% Propriete milieu

D = 202.7e-6 ;  % distance, temps de delai (s) pour d=30cm
C = 1480 ;   % celerite du son dans l'eau (m/s)


%% Propiete emetteur cMUT
 
xie = 0.3;   % amortissement
w0e = 3e6;   % frequence de resonance
numerator = 1;
denominator = [1/w0e^2,2*xie/w0e,1];
syse = tf(numerator,denominator);

%% Propriete recepteur cMUT

xir = 0.3;
w0r = 3e6;
numerator = 1;
denominator = [1/w0r^2,2*xir/w0r,1];
sysr = tf(numerator,denominator);

%% Signal emetteur 

fs = 3e6;  %Frequence du signal
Ncycle = 2;  %Nombre de cycles
% t =0:D/(3*fs):D*1.5;  %D/(3*fs) resolution de l'analyse
% sig_elec = zeros(size(t));
% n = 1;

%generation du signal beaucoup plus rapide que par boucle while
tau = 1/fs;
Tf = D*1.5;
Ts = 1/(50*fs);
[u,t] = gensig("sin",tau,Tf,Ts);
sig_elec=[u(1:round(tau*Ncycle/Ts),1);zeros(length(t)-round(tau*Ncycle/Ts),1)];

% figure;
% plot(t,sig_elec)

% while t(n) <= Ncycle/fs  
%     sig_elec(n) = sin(2*pi*fs*t(n));
%     n=n+1;
% end

% figure;
% plot(t , sig_elec)

sig_acoustic_e = lsim(syse, sig_elec, t); %signal acoustique a l'emission

figure;
plot(t, sig_acoustic_e)

title('signal acoustique émis par le cMUT ')

%% Signal reception t=0

sig_acoustic_r0 = lsim(sysr, sig_acoustic_e, t); %signal acoustique a la reception a t=0

% figure;
% plot(t, sig_acoustic_r0)

% title('signal électrique capté par le cMUT à t=0')
%% Signal reception t=D

sig_acoustic_r = zeros(size(t));

n = 1;
for i =1:length(t)
    if t(i) >= D
        sig_acoustic_r(i) = sig_acoustic_r0(n);
        n = n+1;
    end
end

figure;
plot(t, sig_acoustic_r)
title('signal électrique capté par le cMUT ')

%% Cross correlation

r = xcorr(sig_acoustic_r0,sig_acoustic_r);

figure;
plot(r)
title('correlation croisee')

%a revoir pour determiner le temps de delai
figure;
plot(t,flip(r(1:length(t)))) %marche, mais sans savoir la raison

