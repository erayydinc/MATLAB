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
M_pulse = t / PRI; %% ZAMANDA ÜRETİLEN HER PULSE U PRI YA BÖLDÜK N ADET PULSE ELDE ETTİK.
%n = 2 * max(rd) ./ fd; %%
fd_m = (-N/2:-N/2-1) * fs/N; %% 0 DAN BAŞLAYAN SAMPLİNG


v = 3; R0 = 10;
R1 = R0 + v * t;
tao1 = 2* R1/c;   %% target1 süresi
lambda = f0 /c;

figure (9);

S2_proposed = A1 * rectpuls(fd/ B) .* exp (-1j * 4 * pi / c * (fd + f0) .* R1);
subplot(2,1,1);plot(M_pulse, abs(S2_proposed)); title("S2 proposed");
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

subplot(2,1,2);plot(M_pulse, 20*log10(abs(s2_natural)/max(abs(s2_natural)))); 
title("S2 natural");

%% NUFFT Azimuth Interpolation Parameters
y = 3; % oversampling factor
L = 2; % interpolation length
M = N; % number of azimuth pulses
gamma_M = y * M; % oversampled grid size

% Improved Kaiser-Bessel window function
function phi = kaiser_bessel_window(x, a)
    % Kaiser-Bessel window function using modified Bessel function
    if abs(x) < a
        arg = sqrt(a^2 - x^2);
        phi = besseli(0, L * arg) / besseli(0, L * a);
    else
        phi = 0;
    end
end

% Step 1: Azimuth interpolation (Equation 8)
fprintf('Step 1: Azimuth interpolation...\n');
a = pi * (2 - 1/y) - 0.01; % Kaiser-Bessel parameter

% Initialize interpolated signal
u = zeros(1, gamma_M);

for m = 1:M
    for l = -L:L
        % Non-uniform sampling grid (corresponding to x_m^(n) in equation)
        x_m_n = (1 + fd(m) / f0) * m * PRI;
        
        % Nearest integer to oversampled grid
        mu_m_n = round(y * x_m_n);
        
        % Ensure indices are within bounds
        k_idx = mu_m_n + l;
        if k_idx >= 1 && k_idx <= gamma_M
            % Kaiser-Bessel window argument
            x_arg = y * x_m_n - mu_m_n - l;
            
            % Kaiser-Bessel window function
            phi_val = kaiser_bessel_window(x_arg, a);
            
            % Accumulate interpolated values (Equation 8)
            u(k_idx) = u(k_idx) + (1/sqrt(2*pi)) * phi_val * S2_proposed(m);
        end
    end
end

% Step 2: Compute the azimuth FFT of size γM (Equation 9)
fprintf('Step 2: Computing azimuth FFT...\n');
U = zeros(1, gamma_M);
for m = 1:gamma_M
    for k = -gamma_M/2:(gamma_M/2-1)
        idx = k + gamma_M/2 + 1; % Convert to MATLAB indexing
        if idx >= 1 && idx <= gamma_M
            U(idx) = U(idx) + u(m) * exp(-1j * 2*pi*m*k / gamma_M);
        end
    end
end

% Step 3: Obtain the result via azimuth scaling (Equation 10)
fprintf('Step 3: Azimuth scaling...\n');
S3 = zeros(1, M);
for m = 1:M
    m_scaled = m - M/2; % Scale m to range [-M/2, M/2-1]
    if m_scaled >= -M/2 && m_scaled <= M/2-1
        idx = m_scaled + M/2 + 1; % Convert to MATLAB indexing
        if idx >= 1 && idx <= gamma_M && idx <= length(U)
            % Azimuth scaling factor
            phi_m = kaiser_bessel_window(0, a); % φ(m) at m=0 for scaling
            if phi_m ~= 0
                S3(m) = U(idx) / phi_m;
            end
        end
    end
end

% Display results
figure(10);
subplot(3,1,1);
plot(1:gamma_M, real(u));
title('Step 1: Interpolated signal u(k)');
xlabel('k'); ylabel('Real(u)');

subplot(3,1,2);
plot(1:gamma_M, abs(U));
title('Step 2: Azimuth FFT U(m,n)');
xlabel('m'); ylabel('|U|');

subplot(3,1,3);
plot(1:M, abs(S3));
title('Step 3: Final result S_3(m,n)');
xlabel('m'); ylabel('|S_3|');

fprintf('NUFFT Azimuth Interpolation completed.\n');
fprintf('Original signal length: %d\n', M);
fprintf('Oversampled grid size: %d\n', gamma_M);
fprintf('Interpolation length: %d\n', 2*L+1);