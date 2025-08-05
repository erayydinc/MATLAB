N  = 16;
            % Tp 10 mikro saniyeB  = 1e6;  %% Bandwidth / Tp * B = 1
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
figure(10);
fprintf(num2str(u))
plot(M,real(u));