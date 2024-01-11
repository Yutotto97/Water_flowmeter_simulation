%% Propriete milieu

D = 202.7e-6 ;  % distance, temps de delai (s) pour L=30cm
C = 1480 ;   % celerite du son dans l'eau (m/s)
L=D*C*100 ;  %distance in cm

%% Propiete emetteur cMUT
 
xie = 1.7;   % amortissement
w0e = 3e6;   % frequence de resonance
numerator = 1;
denominator = [1/w0e^2,2*xie/w0e,1];
syse = tf(numerator,denominator);

%% Propriete recepteur cMUT

xir = 1.7;
w0r = 3e6;
numerator = 1;
denominator = [1/w0r^2,2*xir/w0r,1];
sysr = tf(numerator,denominator);

%% Xcorr en fct de Ncycle

fs = 3e6;  %Frequence du signal
Ncycles = 1:10:100;  %Nombre de cycles

%noise_level=0.001; % bruit en % du signal

sigs_e=[];
sigs_r=[];
sigs_xcorr=[];
for Ncycle=Ncycles
    
    tau = 1/fs;
    Tf = D*2;
    Ts = 0.5e-9;
    [u,t] = gensig("sin",tau,Tf,Ts);
    sig_elec=[u(1:round(tau*Ncycle/Ts),1);zeros(length(t)-round(tau*Ncycle/Ts),1)];

    sig_acoustic_e = lsim(syse, sig_elec, t); %signal acoustique a l'emission
    
    sigs_e=[sigs_e,sig_acoustic_e];
    
    sig_acoustic_r0 = lsim(sysr, sig_acoustic_e, t); %signal acoustique a la reception a t=0
    
    sig_acoustic_r = zeros(size(t));

    n = 1;
    for i =1:length(t)
        if t(i) >= D
            sig_acoustic_r(i) = sig_acoustic_r0(n);
            n = n+1;
        end
        sig_acoustic_r(i)=sig_acoustic_r(i)*exp(-0.002*(fs/1e6)^2*L);  %+rand*noise_level;
    end

    sigs_r=[sigs_r,sig_acoustic_r];
    
    r = xcorr(sig_acoustic_r0,sig_acoustic_r,'normalized');
    
    sigs_xcorr=[sigs_xcorr,r];
    
end

figure;
for i=1:length(Ncycles)
    plot(sigs_xcorr(:,i))
    %plot(t,flip(sigs_xcorr(1:length(t),i)))
    hold on
end
legend(split(num2str(Ncycles)));
title('correlation croisee en fct de Ncycles');
xlabel('time (s)');

%Caractérisation du "pic" de la corrélation croisée
%Je relève l'écart de temps entre les 99% du maximum, comme les facteurs de
%qualité

% Qs_N=[];
% for i=1:length(Ncycles)
%     if find(flip(sigs_xcorr(1:length(t),i))==max(sigs_xcorr(1:length(t),i))
%     indices_099=find(flip(sigs_xcorr(1:length(t),i))>=0.999*max(sigs_xcorr(1:length(t),i)));
%     Qs_N=[Qs_N,t(indices_099(end))-t(indices_099(1))];
% end
% figure;
% plot(Ncycles,Qs_N);
% title('Précision en fonction du nombre de cycles');
% xlabel('number of cycles');
% ylabel('time (s)');

%% Xcorr en fct de fs

fss = 1e6:1e6:20e6;  %Frequence du signal

Ncycle = 2;  %Nombre de cycles

noise_level=10; % rapport signal à bruit en dB

sigs_e=[];
sigs_r=[];
sigs_xcorr=[];
errors=[];

for fs=fss
    %def cMUTs emetteur et recepteur
    xie = 1.7;   % amortissement
    w0e = fs;   % frequence de resonance
    numerator = 1;
    denominator = [1/w0e^2,2*xie/w0e,1];
    syse = tf(numerator,denominator);

    xir = 1.7;
    w0r = fs;
    numerator = 1;
    denominator = [1/w0r^2,2*xir/w0r,1];
    sysr = tf(numerator,denominator);
    
%     tau = 1/fs;
%     Tf = D*5;
%     Ts = 1/(50*fss(length(fss)));
%     
    tau = 1/3e6;
    Tf = D*5;
    Ts = 0.5e-9;
    
    [u,t] = gensig("sin",tau,Tf,Ts);
    sig_elec=[u(1:round(tau*Ncycle/Ts),1);zeros(length(t)-round(tau*Ncycle/Ts),1)];

    sig_acoustic_e = lsim(syse, sig_elec, t); %signal acoustique a l'emission
    
    sigs_e=[sigs_e,sig_acoustic_e];
    
    sig_acoustic_r0 = lsim(sysr, sig_acoustic_e, t); %signal acoustique a la reception a t=0
    
    sig_acoustic_r = zeros(size(t));
    
    signal_CC=exp(-0.002*(fs/1e6)^2*L)*(max(sig_acoustic_r0)-min(sig_acoustic_r0));

    noise=6*signal_CC/(10^(noise_level/20)); % amplitude du bruit homogene au signal

    n = 1;
    for i =1:length(t)
        if t(i) >= D
            sig_acoustic_r(i) = sig_acoustic_r0(n);
            n = n+1;
        end
        sig_acoustic_r(i)=sig_acoustic_r(i)*exp(-0.002*(fs/1e6)^2*L)+(rand*2-1)*noise;
    end

    sigs_r=[sigs_r,sig_acoustic_r];
    
    r = xcorr(sig_acoustic_r0,sig_acoustic_r,'normalized');
    
    sigs_xcorr=[sigs_xcorr,r];
    
    answer=t(find(flip(r(1:length(t)))==max(r(1:length(t)))));
    
    errors = [errors,abs(answer-D)/D*100];
    
end

% figure;
% for i=1:length(fss)
%     plot(t(1:300),sigs_e(1:300,i))
%     hold on
% end
% % plot(tx(1:round(length(tx)/2)),sig_philips_Pa(1:round(length(tx)/2)))
% legend(split(num2str(fss)));
% xlabel('time(s)')
% ylabel('Pressure(MPa)')
% title('signal acoustique emise en fonction de la fréquence de résonance');

figure;
for i=1:length(fss)
    %plot(sigs_xcorr(:,i))
    plot(t,flip(sigs_xcorr(1:length(t),i)))
    hold on
end
legend(split(num2str(fss)));
title('correlation croisee en fct de la frequence');

figure;
plot(errors);

% Qs_f=[];
% for i=1:length(fss)
%     indices_099=find(flip(sigs_xcorr(1:length(t),i))>=0.999);
%     Qs_f=[Qs_f,t(indices_099(end))-t(indices_099(1))];
% end
% figure;
% plot(fss,Qs_f);
% title('Précision en fonction de la fréquence');
% xlabel('frequency (Hz)');
% ylabel('time (s)');


% figure;
% for i=1:length(fss)
%     plot(sigs_xcorr(:,i)/max(sigs_xcorr(:,i)))
%     hold on
% end
% legend(split(num2str(fss)));
% title('correlation croisee en fct de frequence (echelle modifiee)');

%% Xcorr en fct du bruit

fs = 3e6;  %Frequence du signal
Ncycle = 2;  %Nombre de cycles

D = 50e-6;  
noise_levels=0:10:100; % rapport signal à bruit en dB

sigs_e=[];
sigs_r=[];
sigs_xcorr=[];
accuracy=[];

for noise_level=noise_levels
    tau = 1/fs;
    Tf = D*1.5;
    Ts = 0.5e-9;
    [u,t] = gensig("sin",tau,Tf,Ts);
    sig_elec=[u(1:round(tau*Ncycle/Ts),1);zeros(length(t)-round(tau*Ncycle/Ts),1)];

    sig_acoustic_e = lsim(syse, sig_elec, t); %signal acoustique a l'emission
    
    sigs_e=[sigs_e,sig_acoustic_e];
    
    sig_acoustic_r0 = lsim(sysr, sig_acoustic_e, t); %signal acoustique a la reception a t=0
    
    sig_acoustic_r = zeros(size(t));

    n = 1;
    
    signal_CC=exp(-0.002*(fs/1e6)^2*L)*(max(sig_acoustic_r0)-min(sig_acoustic_r0))
    
    noise=6*signal_CC/(10^(noise_level/20)) % amplitude du bruit homogene au signal

    for i =1:length(t)
        if t(i) >= D
            sig_acoustic_r(i) = sig_acoustic_r0(n);
            n = n+1;
        end
        sig_acoustic_r(i)=sig_acoustic_r(i)*exp(-0.002*(fs/1e6)^2*L)+(rand*2-1)*noise;
    end
    
    sigs_r=[sigs_r,sig_acoustic_r];
    
    r = xcorr(sig_acoustic_r0,sig_acoustic_r,'normalized');
    
    answer=t(find(flip(r(1:length(t)))==max(r(1:length(t)))));
    
    accuracy=[accuracy,answer];
    
    sigs_xcorr=[sigs_xcorr,r];
    
end

figure;
subplot(2,1,1)
plot(t,sig_acoustic_r0)
hold on;
plot(t,sig_acoustic_r)
xline(D,'--')
title('emitted and received signal for f=3MHz, xi=1.7, SNR=0dB');
legend('emitted signal','received signal',strcat('Tof= ',num2str(D)));
xlabel('time (s)');

subplot(2,1,2)
plot(t,flip(sigs_xcorr(1:length(t),end)))
xline(t(find(flip(r(1:length(t)))==max(r(1:length(t))))),'--');
title('cross correlation for f=3MHz, xi=1.7, SNR=0dB');
legend('cross correlation', strcat('Tof=',num2str(t(find(flip(r(1:length(t)))==max(r(1:length(t))))))))
xlabel('time (s)');

figure;
for i=1:length(noise_levels)
    %plot(sigs_xcorr(:,i))
    plot(t,flip(sigs_xcorr(1:length(t),i)))
    hold on
end
legend(split(num2str(noise_levels)));
title('cross correlation as a function of SNR (dB)');
xlabel('time (s)');

%Caractérisation du "pic" de la corrélation croisée
%Je relève l'écart de temps entre les 99% du maximum, comme les facteurs de
%qualité

% Qs_SNR=[];
% for j=1:length(noise_levels)
%     answer=t(find(flip(sigs_xcorr(1:length(t),j))==max(sigs_xcorr(1:length(t),j))));
%     if answer>=D*0.999 & answer<=D*1.001
%         indices_099=find(flip(sigs_xcorr(1:length(t),j))>=0.999*max(sigs_xcorr(1:length(t),j)));
%         Qs_SNR=[Qs_SNR,t(indices_099(end))-t(indices_099(1))];
%     else
%         Qs_SNR=[Qs_SNR,0];
%     end
% end
% figure;
% plot(noise_levels,Qs_SNR);
% title("Précision en fonction du SNR");
% xlabel('SNR (dB)');
% ylabel('time (s)');

figure;
plot(noise_levels,accuracy);
yline(D);
title("temps de transit calculé en fonction du SNR");
xlabel('bruit SNR (dB)');
ylabel('temps (s)');

%% Xcorr par répétition

fs = 3e6;  %Frequence du signal
Ncycle = 2;  %Nombre de cycles

noise_level=40; % rapport signal à bruit en dB

repetitions=1:5:101; %obtention du résultat en fonction du nombre de répétitions

sigs_e=[];
sigs_r=[];
sigs_xcorr=[];
accuracy=[];

for repetition=repetitions
    tau = 1/fs;
    Tf = D*2;
    Ts = 0.5e-9;
    [u,t] = gensig("sin",tau,Tf,Ts);
    sig_elec=[u(1:round(tau*Ncycle/Ts),1);zeros(length(t)-round(tau*Ncycle/Ts),1)];

    sig_acoustic_e = lsim(syse, sig_elec, t); %signal acoustique a l'emission
    
    sigs_e=[sigs_e,sig_acoustic_e];
    
    sig_acoustic_r0 = lsim(sysr, sig_acoustic_e, t); %signal acoustique a la reception a t=0
    
    sig_acoustic_rs = zeros(size(t));
    sig_acoustic_r = zeros(size(t));
    
    answers=0.;
    for rep=1:repetition
        n = 1;
        noise=6/(sqrt(2)*10^(noise_level/20)); % conversion du bruit homogene au signal
        for i =1:length(t)
            if t(i) >= D
                sig_acoustic_r(i) = sig_acoustic_r0(n);
                n = n+1;
            end
            sig_acoustic_r(i)=sig_acoustic_r(i)*exp(-0.002*(fs/1e6)^2*L)+(rand*2-1)*noise;
        end
        sig_acoustic_rs=sig_acoustic_rs+sig_acoustic_r;
        r = xcorr(sig_acoustic_r0,sig_acoustic_r,'normalized');
        answers=answers+t(find(flip(r(1:length(t)))==max(r(1:length(t)))));
    end
    
    ToF=answers/repetition;
    
    accuracy=[accuracy,abs(ToF-D)/D*100];
    
    sig_acoustic_rs=sig_acoustic_rs/repetition;
    
    sigs_r=[sigs_r,sig_acoustic_rs];
    
    r = xcorr(sig_acoustic_r0,sig_acoustic_r,'normalized');
    
    sigs_xcorr=[sigs_xcorr,r];
    
end

% figure;
% for i=1:length(repetitions)
%     %plot(sigs_xcorr(:,i))
%     plot(t,flip(sigs_xcorr(1:length(t),i)))
%     hold on
% end
% legend(split(num2str(repetitions)));
% xline(D)
% title('correlation croisee en fct du nombre de répétition');
% xlabel('time (s)');
% 
% %Caractérisation du "pic" de la corrélation croisée
% %Je relève l'écart de temps entre les 99% du maximum, comme les facteurs de
% %qualité
% 
% Qs_rep=[];
% for j=1:length(repetitions)
%     answer=t(find(flip(sigs_xcorr(1:length(t),j))==max(sigs_xcorr(1:length(t),j))));
%     if answer>=D*0.999 & answer<=D*1.001
%         indices_0999=find(flip(sigs_xcorr(1:length(t),j))>=0.999*max(sigs_xcorr(1:length(t),j)));
%         Qs_rep=[Qs_rep,t(indices_0999(end))-t(indices_0999(1))];
%     else
%         Qs_rep=[Qs_rep,0];
%     end
% end
% figure;
% plot(repetitions,Qs_rep);
% title("Précision en fonction du nombre de répétition");
% xlabel('répétition');
% ylabel('time (s)');

figure;
plot(repetitions,accuracy);
title("Erreur en fonction du nombre de répétition (SNR=40dB)");
xlabel('répétition');
ylabel('Erreur en %');

%% Xcorr moyenné en fonction de la frequence

Ncycle = 2;  %Nombre de cycles

noise_level=40; % rapport signal à bruit en dB

repetition=50; %obtention du résultat en fonction du nombre de répétitions

fss = 1e6:1e6:20e6;  %Frequence du signal

sigs_e=[];
sigs_r=[];
sigs_xcorr=[];
accuracy=[];

for fs=fss
    %def cMUTs emetteur et recepteur
    xie = 1.7;   % amortissement
    w0e = fs;   % frequence de resonance
    numerator = 1;
    denominator = [1/w0e^2,2*xie/w0e,1];
    syse = tf(numerator,denominator);

    xir = 1.7;
    w0r = fs;
    numerator = 1;
    denominator = [1/w0r^2,2*xir/w0r,1];
    sysr = tf(numerator,denominator);
    
    tau = 1/fs;
    Tf = D*2;
    Ts = 0.5e-9;
    [u,t] = gensig("sin",tau,Tf,Ts);
    sig_elec=[u(1:round(tau*Ncycle/Ts),1);zeros(length(t)-round(tau*Ncycle/Ts),1)];

    sig_acoustic_e = lsim(syse, sig_elec, t); %signal acoustique a l'emission
    
    sigs_e=[sigs_e,sig_acoustic_e];
    
    sig_acoustic_r0 = lsim(sysr, sig_acoustic_e, t); %signal acoustique a la reception a t=0
    
    sig_acoustic_rs = zeros(size(t));
    sig_acoustic_r = zeros(size(t));
    answers=0.;
    for rep=1:repetition
        n = 1;
        noise=6/(sqrt(2)*10^(noise_level/20)); % conversion du bruit homogene au signal
        for i =1:length(t)
            if t(i) >= D
                sig_acoustic_r(i) = sig_acoustic_r0(n);
                n = n+1;
            end
            sig_acoustic_r(i)=sig_acoustic_r(i)*exp(-0.002*(fs/1e6)^2*L)+(rand*2-1)*noise;
        end
        sig_acoustic_rs=sig_acoustic_rs+sig_acoustic_r;
        r = xcorr(sig_acoustic_r0,sig_acoustic_r,'normalized');
        answers=answers+t(find(flip(r(1:length(t)))==max(r(1:length(t)))));
    end
    
    ToF=answers/repetition;
    
    accuracy=[accuracy,abs(ToF-D)/D*100];
    
    sig_acoustic_rs=sig_acoustic_rs/repetition;
    
    sigs_r=[sigs_r,sig_acoustic_rs];
    
    r = xcorr(sig_acoustic_r0,sig_acoustic_r,'normalized');
    
    sigs_xcorr=[sigs_xcorr,r];
    
end

figure;
for i=1:length(fss)
    %plot(sigs_xcorr(:,i))
    plot(t,flip(sigs_xcorr(1:length(t),i)))
    hold on
end
legend(split(num2str(fss)));
xline(D)
title('correlation croisee en fct du nombre de répétition');
xlabel('time (s)');

%Caractérisation du "pic" de la corrélation croisée
%Je relève l'écart de temps entre les 99% du maximum, comme les facteurs de
%qualité

Qs_rep_fs=[];
for j=1:length(fss)
    answer=t(find(flip(sigs_xcorr(1:length(t),j))==max(sigs_xcorr(1:length(t),j))));
    if answer>=D*0.999 & answer<=D*1.001
        indices_0999=find(flip(sigs_xcorr(1:length(t),j))>=0.999*max(sigs_xcorr(1:length(t),j)));
        Qs_rep_fs=[Qs_rep_fs,t(indices_0999(end))-t(indices_0999(1))];
    else
        Qs_rep_fs=[Qs_rep_fs,0];
    end
end

figure;
plot(fss,Qs_rep_fs);
title("Précision en fonction de la fréquence (SNR=40dB)");
xlabel('fréquence (MHz)');
ylabel('time (s)');

figure;
plot(fss,accuracy);
title("Erreur en fonction de la fréquence (SNR=40dB)");
xlabel('fréquence (MHz)');
ylabel('Erreur en %');

%% Limite de distance en focntion de la fréquence

D = 202.7e-6 ;  % distance, temps de delai (s) pour L=30cm
C = 1480 ;   % celerite du son dans l'eau (m/s)
L=D*C*100 ;  %distance in cm

Ncycle = 2;  %Nombre de cycles

noise_level=10; % rapport signal à bruit en dB

%signal CC et amplitude CC du bruit pour 10dB, 3MHz, 30cm :
%signal_CC = 0.0118;
noise_CC = 0.0159;

repetition=100; %obtention du résultat en fonction du nombre de répétitions
% 
% fss = 1e6:1e6:20e6;  %Frequence du signal
% 
% Ls=1:1:100;  %distance en cm

error_limit=1;  %erreur maximale tolérée en %

fss = 1e6:1e6:20e6;  %Frequence du signal
Ls=1:1:100;  %distance en cm

limits=[];

signal_CCs=[];

figure;

h=animatedline;
axis([fss(1) fss(end) Ls(1) Ls(end)]);

for fs=fss
    n=1;
    accuracy=0.;
    %def cMUTs emetteur et recepteur
    xie = 1.7;   % amortissement
    w0e = fs;   % frequence de resonance
    numerator = 1;
    denominator = [1/w0e^2,2*xie/w0e,1];
    syse = tf(numerator,denominator);

    xir = 1.7;
    w0r = fs;
    numerator = 1;
    denominator = [1/w0r^2,2*xir/w0r,1];
    sysr = tf(numerator,denominator);
    
    while n<length(Ls) & accuracy<error_limit
        D=Ls(n)/(C*100);
        tau = 1/fs;
        Tf = D*2;
        Ts = 0.5e-9;
        [u,t] = gensig("sin",tau,Tf,Ts);
        sig_elec=[u(1:round(tau*Ncycle/Ts),1);zeros(length(t)-round(tau*Ncycle/Ts),1)];

        sig_acoustic_e = lsim(syse, sig_elec, t); %signal acoustique a l'emission

        sig_acoustic_r0 = lsim(sysr, sig_acoustic_e, t); %signal acoustique a la reception a t=0

        sig_acoustic_rs = zeros(size(t));
        sig_acoustic_r = zeros(size(t));
        
        signal_CC=exp(-0.002*(fs/1e6)^2*L)*(max(sig_acoustic_r0)-min(sig_acoustic_r0));
%     
%         noise_CC=6*signal_CC/(sqrt(2)*10^(noise_level/20));
    
        answers=0.;
        
        for rep=1:repetition   
            j = 1;
            for i=1:length(t)
                if t(i) >= D
                    sig_acoustic_r(i) = sig_acoustic_r0(j);
                    j = j+1;
                end
                sig_acoustic_r(i)=sig_acoustic_r(i)*exp(-0.002*(fs/1e6)^2*Ls(n))+(rand*2-1)*noise_CC/2;
            end
            sig_acoustic_rs=sig_acoustic_rs+sig_acoustic_r;
            r = xcorr(sig_acoustic_r0,sig_acoustic_r,'normalized');
            answers=answers+t(find(flip(r(1:length(t)))==max(r(1:length(t)))));
        end
        ToF=answers/repetition;
        accuracy=abs(ToF-D)/D*100;
        if accuracy<error_limit
            n=n+1;
        end
    end
    
    limits=[limits, Ls(n)];
    
    signal_CCs=[signal_CCs,signal_CC];
    addpoints(h,fs,Ls(n));
    
    drawnow
end

figure;
plot(fss,limits);
title("Limite de distance en fonction de la fréquence (SNR=10dB, Erreur<1%)");
xlabel('fréquence (MHz)');
ylabel('Distance');

limit_SNRs=zeros(length(signal_CCs),1);
n=1;
for scc=signal_CCs
    limit_SNRs(n)=20*log(scc/(sqrt(2)*noise_CC/6));
    n=n+1;
end
figure;
plot(fss,limit_SNRs);
title("tolerated SNR as a function of frequency")
xlabel('frequency (MHz)');
ylabel('SNR (dB)');

