% Part 2.5: use DFT to find frequency burried in noise and amplitudes of the peak frequencies using Fourier transform
% signal:

S = 1.8 + 0.7*sin(2*pi*70*t) + sin(2*pi*130*t);

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
