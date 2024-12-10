% Part 1A:
% y1(n) = x(n) - x(n-6)
Fs = 360;

a = 1;
b1 = [1 0 0 0 0 0 -1];

[h1, w1] = freqz(b1, a);
f1 = w1/pi * (Fs/2);

figure;
subplot(311);
stem(f1, abs(h1));
title('Magnitude of the frequency response of y1(n) filter');

Ts = 1/Fs;
L = 1500; % Length of the signal
t = (0:L-1)/Fs; % Time vector

x = cos(120*pi*t) + cos (240*pi*t) + 3 + cos(60*pi*t);
y1 = filter(b1, a, x);

% UNFILTERED
Y_uf1 = fft(x, L); % FFT for unflitered signal
f_uf_fft1 = Fs/L * (-L/2:L/2-1); % frequency vector 

subplot(312);
plot(f_uf_fft1, abs(fftshift(Y_uf1)), 'LineWidth', 3);
title('FFT Spectrum of the Unfiltered Signal');
xlabel('Frequency (Hz)');
ylabel('|FFT(Y)|');

% FILTERED
Yk1 = fft(y1, L); % FFT of the filtered signal
f_fft1 = Fs/L * (-L/2:L/2-1); % Frequency vector

subplot(313);
plot(f_fft1, abs(fftshift(Yk1)), 'LineWidth', 3);
title('FFT Spectrum of the Filtered Signal');
xlabel('Frequency (Hz)');
ylabel('|FFT(Y)|');

% plot sinusoidal signal
figure;
plot(t, x, t, y1);
title('Original and filtered signals');
legend('Original', 'Filtered');
xlabel('Time (s)');
ylabel('Amplitude');

% Part 1B:
% y2(n) = x(n) + x(n-6)
a = 1;
b2 = [1 0 0 0 0 0 1];

[h2, w2] = freqz(b2, a);
f2 = w2/pi * (Fs/2);

figure;
subplot(311);
stem(f2, abs(h2));
title('Magnitude of the frequency response of y2(n) filter');

y2 = filter(b2, a, x);

% UNFILTERED
Y_uf2 = fft(x, L); % FFT for unflitered signal
f_uf_fft2 = Fs/L * (-L/2:L/2-1); % frequency vector 

subplot(312);
plot(f_uf_fft2, abs(fftshift(Y_uf2)), 'LineWidth', 3);
title('FFT Spectrum of the Unfiltered Signal');
xlabel('Frequency (Hz)');
ylabel('|FFT(Y)|');

% FILTERED
Yk2 = fft(y2, L); % FFT of the filtered signal
f_fft2 = Fs/L * (-L/2:L/2-1); % Frequency vector

subplot(313);
plot(f_fft2, abs(fftshift(Yk2)), 'LineWidth', 3);
title('FFT Spectrum of the Filtered Signal');
xlabel('Frequency (Hz)');
ylabel('|FFT(Y)|');

% Part 2: Discrete Fourier Transform
xn = [0, 1, 2, 3];
N = length(xn);

n = [0: 1 :N-1]; % row vector for n
k= [0: 1: N-1]; % row vector for k

WN = exp(-j*2*pi/N); % Wn factor
nk = n'* k; % Creates an N by N matrix of nk values
WNnk = WN .^ nk ; % row vector for DEFT coefficients

Xk = xn * WNnk; % compute DFT coefficient
Xk_fft = fft(xn);

disp('DFT:');
disp(Xk_fft);

% Part 2.2: plot magnitude using stem command
figure;
subplot(211);
stem(n, abs(Xk));
title('Magnitude of Discrete Fourier Transform of x(n)');
xlabel('n');
ylabel('Magnitude');

% Part 2.3: Repeat 2 using the FFT algorithm (refer to the appendix)
subplot(212);
stem(0:N-1, abs(Xk_fft))
title('Magnitude of Discrete Fourier Transform of x(n) using FFT');
xlabel('n');
ylabel('Magnitude');

% Part 2.4: Using X(k), calculate the inverse discrete Fourier transform and obtain x(n)
X_ifft = ifft(Xk_fft);
disp('Inverse DDT transform:');
disp(X_ifft);

% Part 2.5: use DFT to find frequency burried in noise and amplitudes of the peak frequencies using Fourier transform
% signal:

S = 0.8 + 0.7*sin(2*pi*50*t) + sin(2*pi*120*t);

Fs2 = 1000; % Sampling frequency

T = 1/Fs; % Sampling period
L = 1500; % Length of signal
t = (0:L-1)*T; % Time vector

X_buried = S + 2*randn(size(t));

% Corrupt the signal with zero-mean random noise with a variance of 4.
figure;
plot(1000*t,X_buried)
title("Signal Corrupted with Zero-Mean Random Noise")
xlabel("t (milliseconds)")
ylabel("X(t)")

% Compute the Fourier transform of the signal.
Y_buried = fft(X_buried);
plot(Fs/L*(-L/2:L/2-1),abs(fftshift(Y_buried)),"LineWidth",3)
title("fft Spectrum in the Positive and Negative Frequencies")
xlabel("f(Hz)")
ylabel("|fft(X)|")

% find the amplitudes of the peak frequencies using Fourier transform
P2 = abs(Y_buried/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

% Plot the single-sided amplitude spectrum P1
f = Fs/L*(0:(L/2));
figure;
plot(f,P1,"LineWidth",3) 
title("Single-Sided Amplitude Spectrum of X(t)")
xlabel("f (Hz)")
ylabel("|P1(f)|")

% Fourier transform of the original, uncorrupted signal and retrieve the exact amplitudes
Y_signal = fft(S);
P2 = abs(Y_signal/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

figure;
plot(f,P1,"LineWidth",3) 
title("Single-Sided Amplitude Spectrum of S(t)")
xlabel("f (Hz)")
ylabel("|P1(f)|")

% print out peak amplitude
[pks, locs] = findpeaks(P1, f, 'MinPeakHeight', 0.1); % Adjust MinPeakHeight as needed
disp('Peak Frequencies (Hz):');
disp(locs);
disp('Amplitudes of Peak Frequencies:');
disp(pks);
