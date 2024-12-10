% signal frequency and period of the signal 1
F1 = 100;
T1 = 1/F1;

Fs = 50; % sampling frequency
Ts = 1/Fs; % sampling period 

t = 0: T1/100 : 3*T1; % signal time axis

xc1 = 10*cos(200*pi*t); % signal 1
xc2 = 10*cos(100*pi*t); % signal 2

tx = 0: Ts :3*Ts; % sampling time axis
xs = 10*cos(200*pi*tx);

% plot trong cung mot axis
subplot(211)
plot(t*2, xc1, 'LineWidth', 2);
hold on
plot(t*2, xc2, 'LineWidth', 2);
stem(tx, xs, "LineWidth", 2);

Fs2 = 500;
Ts2 = 1/Fs2;

tx2 = 0: Ts2: 3*T1;
xs2 = 10*cos(200*pi*tx2); % sampling signal

subplot(212)
plot(t, xc1, 'LineWidth', 2);
hold on
plot(t, xc2, 'LineWidth', 2);
stem(tx2, xs2, "LineWidth", 2);
title('Sampling with 50Hz and 500Hz');
hold off

% Part 2
% 2.1. Impulse
% FIR
n = -10:1:39;
a = 1; % idk what this shit is
b = [0.2, 0.2, 0.2, 0.2, 0.2]; % coefficient
x = [zeros(1,10), 1, zeros(1,39)];
y = filter(b, a, x);

% all impulse plots in 1 figure:
figure(3);
subplot(221);
stem(y);
title('FIR for impulse');

subplot(222);
stem(n, y);
title('FIR index n for impulse')

% IIR
a2 = [1 -1];
b2 = [0.2 0 0 0 0 -0.2];
y2 = filter(b2, a2, x);

subplot(223);
stem(y2);
title('IIR for impulse');

subplot(224);
stem(n, y2);
title('IIR with n for impulse')

% 2.2. Step 
% FIR 
x3 = [zeros(1,10), ones(1,40)];
y3 = filter(b, a, x3);

figure(4);
subplot(221);
stem(y3);
title('FIR for step')

subplot(222);
stem(n, y3);
title('FIR index n for step')

% IIR
y4 = filter(b2, a2, x3);

subplot(223);
stem(y4);
title('IIR for step');

subplot(224);
stem(n, y4);
title('IIR index n for step');

%convolution a(i) - FIR

w1 = conv(x, y);
w1 = w1(1:length(n));
figure(5);
stem(n, w1, 'filled'); 
title("a(i)");

% convolution a(ii) - unit system impulse responses h(n)

w2 = conv(x, y2);
w2 = w2(1:length(n));
figure(6);
stem(n, w2, "filled")
title("a(ii)")

% Part B:
% Frequency
% using on FIR
figure(7);
subplot(211);
[h,w] = freqz(b, a, 'whole', 2001);
plot(w/pi,20*log10(abs(h)));

% IIR
display(a2);
subplot(212);
[h,w] = freqz(b2, a2, 'whole', 2001);
plot(w/pi,20*log10(abs(h)));
