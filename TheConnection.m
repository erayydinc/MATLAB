% S2 yi oluşturma deneyi//
%% S1 and s2 From the Real Srecieved
N  = 4000;
Tp = 10e-6; %% Tp 10 mikro saniye
B  = 10e6;  %% Bandwidth / Tp * B = 1
f0 = 10e9;  %% Center Frequency -- carrier frequency
fs = 100e6; %% Fsampling - frekasn domainde baktığımız aralık

c  = 3e8;
t  = (0:N-1)/fs;  %% N/fs normalde pulse width eder burada 4 pulse width seçilmiş
rd = c/2 .* t; %% zamana bağlı alınan yol
fd = (-N/2:N/2-1)*fs/N;%% -fs/2 to fs/2
kr = B/Tp; %% chirp rate
A1 = 1;

PRI = 20e-6;
m = 2 * t / PRI; %% pulse sayısı (N adet)
n = 2 * max(rd) ./ fd; %%

v = 100; R0 = 3000;
R1 = R0 + v * m;
tao1 = 2* R1/c;   %% target1 süresi
lambda = f0 /c;




sr  = rectpuls(t-Tp/2-tao1,Tp).*exp(1i*pi*kr*(t-Tp/2-tao1).^2).*exp(-1i * 4 * pi * R1 / lambda);
sf  = rectpuls(t-Tp/2,Tp).*exp(1i*pi*kr*(t-Tp/2).^2);
sf_h =rectpuls(t-Tp/2,Tp).*exp(1i*pi*kr*(t-Tp/2).^2).*[hamming(1001)',zeros(1,(N-1001))];

sr_fft = fft(sr,N); %%
sf_fft = fft(sf,N); %%
sf_hfft = fft(sf_h,N);
s_out = ifft(sr_fft.*conj(sf_hfft)); %% frekans domaininde çarpım.

s1_natural = s_out;
s2_natural = fftshift(fft(s1_natural));

figure(1);

subplot(321);plot(rd/1e3,real(sr)); %% real of s recieved
xlabel('Range/km');ylabel('Amplitude');title('Signal');
subplot(322);plot(fd/1e6,abs(fftshift(sr_fft))); %% Fourier of s recieved
xlabel('Frequency/MHz');ylabel('Amplitude');title("SIGNAL RECIEVED");

subplot(3,2,[5 6]);
plot(rd/1e3,20*log10(abs(s_out)/max(abs(s_out))));
xlabel('Range/km');ylabel('Amplitude/dB');title("SIGNAL COMPRESSED/S1 NATURAL");
axis ([0 6 -60 1]);

figure(2);
subplot(3,2,[1 2]);
plot(fd,20*log10(abs(s2_natural)/max(abs(s2_natural))));
xlabel('Frequency');ylabel('Amplitude/dB');title("S2 NATURAL");
subplot(3,2,[3 4]);
plot(rd/1e3,20*log10(abs(s2_natural)/max(abs(s2_natural))));
xlabel('Range/km');ylabel('Amplitude/dB');title("S2 NATURAL");

%% S1 and S2 From the Proposed Formulas

figure(3);
S1_proposed = sinc(B .* (t-tao1)) .* exp(-1j * 4 * pi * R1 / lambda);

subplot(3,1,1); plot(t*1e6,abs(S1_proposed));
xlabel("time in micro second");
ylabel("abs S1_Proposed");
subplot(3,1,2);
plot(rd/1e3,20*log10(abs(S1_proposed)/max(abs(S1_proposed))));
xlabel('Range/km');ylabel('Amplitude/dB');title("S1 Proposed FORMULA");
axis ([0 6 -60 1]);

S2 = fftshift(fft(S1_proposed));
figure(4);
subplot(3,1,1); plot(fd,abs(S2));
xlabel("frequency");
ylabel("abs S2_Present");
subplot(3,1,2);
plot(rd/1e3,20*log10(abs(S2)/max(abs(S2))));
xlabel('Range/km');ylabel('Amplitude/dB');title("S2 CONVERTED FORMULA");
axis ([0 6 -60 1]);

S2_proposed = A1 * rectpuls(fd / B) .* exp (-1j * 4 * pi / c * (fd + f0) .* R1);
figure(5);
subplot(3,1,1); plot(fd,abs(S2_proposed));
xlabel("frequency");
ylabel("abs S2_proposed");
subplot(3,1,2);
plot(rd/1e3,20*log10(abs(S2_proposed)/max(abs(S2_proposed))));
xlabel('Range/km');ylabel('Amplitude/dB');title("S2_proposed");
axis ([0 6 -60 1]);

%% Şimdi 2 yol var 1. S2_NATURAL 2. S2_PROPOSED KULLANARAK YAZMAK.

%% S2_NATURAL(m,n)

% MATLAB Ayrık yapma kodu kodu OLMADI

sr  = rectpuls(m-Tp/2-tao1,Tp).*exp(1i*pi*kr*(m-Tp/2-tao1).^2).*exp(-1i * 4 * pi * R1 / lambda);
sf  = rectpuls(m-Tp/2,Tp).*exp(1i*pi*kr*(m-Tp/2).^2);
sf_h =rectpuls(m-Tp/2,Tp).*exp(1i*pi*kr*(m-Tp/2).^2).*[hamming(999)',zeros(1,(N-999))];

sr_fft = fft(sr,N); %%
sf_fft = fft(sf,N); %%
sf_hfft = fft(sf_h,N);
s_out = ifft(sr_fft.*conj(sf_hfft)); %% frekans domaininde çarpım.

s1_natural = s_out;
s2_natural_discrete = fftshift(fft(s1_natural));  %% s2 yi tm ile yazıcaz her bir pulse u göstericek




%% GEMİNİ KONDUNDA KALDIM


figure(6);
subplot(3,2,[1 2]);
plot(m,20*log10(abs(s2_natural)/max(abs(s2_natural))));
xlabel('Pulse Number');ylabel('Amplitude/dB');title("S2 NATURAL");
subplot(3,2,[3 4]);
plot(n,20*log10(abs(s2_natural)/max(abs(s2_natural))));
xlabel('Range Bin');ylabel('Amplitude/dB');title("S2 NATURAL");




%% S3 sinyali plotting
N  = 4000;
Tp = 10e-6; %% Tp 10 mikro saniye
B  = 10e6;  %% Bandwidth / Tp * B = 1
f0 = 10e9;  %% Center Frequency -- carrier frequency
fs = 100e6; %% Fsampling - frekasn domainde baktığımız aralık

c  = 3e8;
t  = (0:N-1)/fs;  %% N/fs normalde pulse width eder burada 4 pulse width seçilmiş
rd = c/2 .* t; %% zamana bağlı alınan yol
fd = (-N/2:N/2-1)*fs/N;%% -fs/2 to fs/2
kr = B/Tp; %% chirp rate
A1 = 1;
R0 = 100;
A3 = 1;
fd_doppler = max(fd) / 2;
r = rectpuls(fd, B);
e = exp(-1i * 4 * pi / c * (fd + f0) * R0);

S3 = A3 * r .* e ;
figure (7);
plot(fd, real(S3));

%%
clear all;

N  = 4000;
Tp = 10e-6; %% Tp 10 mikro saniye
B  = 10e6;  %% Bandwidth / Tp * B = 1
f0 = 10e9;  %% Center Frequency -- carrier frequency
fs = 100e6; %% Fsampling - frekasn domainde baktığımız aralık

c  = 3e8;
t  = (0:N-1)/fs;  %% N/fs normalde pulse width eder burada 4 pulse width seçilmiş
rd = c/2 .* t; %% zamana bağlı alınan yol
fd = (-N/2:N/2-1)*fs/N;%% -fs/2 to fs/2
kr = B/Tp; %% chirp rate
A1 = 1;

PRI = 1/100 *1e-6;
m = t / PRI ;


%% ZAMANDA ÜRETİLEN HER PULSE U PRI YA BÖLDÜK N ADET PULSE ELDE ETTİK.
%n = 2 * max(rd) ./ fd; %%

v = 100; R0 = 3000;
R1 = R0 + v * t;
tao1 = 2* R1/c;   %% target1 süresi
lambda = f0 /c;

figure (8);

S2_proposed = A1 * rectpuls(fd / B) .* exp (-1j * 4 * pi / c * (fd + f0) .* R1);
subplot(2,1,1);plot(m, 20*log10(abs(S2_proposed)/max(abs(S2_proposed)))); title("S2 proposed");

sr  = rectpuls(t-Tp/2-tao1,Tp).*exp(1i*pi*kr*(t-Tp/2-tao1).^2).*exp(-1i * 4 * pi * R1 / lambda);
sf  = rectpuls(t-Tp/2,Tp).*exp(1i*pi*kr*(t-Tp/2).^2);
sf_h =rectpuls(t-Tp/2,Tp).*exp(1i*pi*kr*(t-Tp/2).^2).*[hamming(1001)',zeros(1,(N-1001))];

sr_fft = fft(sr,N); %%
sf_fft = fft(sf,N); %%
sf_hfft = fft(sf_h,N);
s_out = ifft(sr_fft.*conj(sf_hfft)); %% frekans domaininde çarpım.

s1_natural = s_out;
s2_natural = fftshift(fft(s1_natural));

subplot(2,1,2);plot(m, 20*log10(abs(s2_natural)/max(abs(s2_natural)))); 
title("S2 natural");

%% Lower the pulse number, replot the grapph
% make the kaiser window % give value to variables (oversampling factor
% etc.)
% then build the interpolation

clear all;

N  = 16;
            % Tp 10 mikro saniye
B  = 1e6;  %% Bandwidth / Tp * B = 1
f0 = 1e9;  %% Center Frequency -- carrier frequency
fs = 10e6; %% Fsampling - frekasn domainde baktığımız aralık

c  = 3e8;
t  = (0:N-1)/fs;  %% N/fs normalde pulse width eder burada 4 pulse width seçilmiş
rd = c/2 .* t; %% zamana bağlı alınan yol
fd = (-N/2:N/2-1)*fs/N;%% -fs/2 to fs/2

A1 = 1;

PRI = 1/fs;
Tp = PRI;
kr = B/Tp; %% chirp rate
M = t / PRI; %% ZAMANDA ÜRETİLEN HER PULSE U PRI YA BÖLDÜK N ADET PULSE ELDE ETTİK.
%n = 2 * max(rd) ./ fd; %%
fd_m = (-N/2:-N/2-1) * fs/N; %% 0 DAN BAŞLAYAN SAMPLİNG


v = 3; R0 = 10;
R1 = R0 + v * t;
tao1 = 2* R1/c;   %% target1 süresi
lambda = f0 /c;

figure (9);

S2_proposed = A1 * rectpuls(fd/ B) .* exp (-1j * 4 * pi / c * (fd + f0) .* R1);
subplot(2,1,1);plot(M, abs(S2_proposed)); title("S2 proposed");
axis ([ 1 64 -2 2]);

sr  = rectpuls(t-Tp/2-tao1,Tp).*exp(1i*pi*kr*(t-Tp/2-tao1).^2).*exp(-1i * 4 * pi * R1 / lambda);
sf  = rectpuls(t-Tp/2,Tp).*exp(1i*pi*kr*(t-Tp/2).^2);
sf_h =rectpuls(t-Tp/2,Tp).*exp(1i*pi*kr*(t-Tp/2).^2).*[hamming(N/4)',zeros(1,(N-N/4))];

sr_fft = fft(sr,N); %%
sf_fft = fft(sf,N); %%
sf_hfft = fft(sf_h,N);
s_out = ifft(sr_fft.*conj(sf_hfft)); %% frekans domaininde çarpım.

s1_natural = s_out;
s2_natural = fftshift(fft(s1_natural));

subplot(2,1,2);plot(M, 20*log10(abs(s2_natural)/max(abs(s2_natural)))); 
title("S2 natural");

%% kaiser window
y = 3; %oversampling factor
L = 2; % interpolation lenght
M = t / PRI;

function y_1 = window(x)
    
    syms k;
    seri = ((x/2) .^ k/ factorial(k)) .^ 2;
    toplam =  symsum(seri, k, 0, inf);
    y_1 = toplam;
end
 % Non uniform sampling grid
a = pi * (2 - 1/y ) - 0.01;

u = 0;


for m = 1:N
    for l = -L:L 
       x_m = (1 + fd(m) / f0) .* m * PRI;
       mu = round(y*x_m);
       x = L * sqrt(a^2 -(y*x_m - mu - l)^ 2 );
       kaiserBessel = window(x);
    
       if abs(x) < a
            teta = kaiserBessel;
        else 
            teta = 0;
        end    
       u = u + 1 / sqrt (2 * pi) * teta * S2_proposed(m);
    end
end

fprintf(num2str(u))
plot(M,real(u));