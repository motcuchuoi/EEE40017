t = linspace(-pi,pi,100);
rng default % initialize random number generator
x = sin(t) + 0.25*rand(size(t));

figure;
plot(t, x);
title('Sinusoidal signal corrupted by random noise implementation');

% filter: y(n)=1/windowSize(x(n)+x(n−1)+...+x(n−(windowSize−1))).
% use different window sizes 5,7 and 12

a = 1;
b1 = [0.2 0.2 0.2 0.2 0.2]; % size 5

b2 = [1/7 1/7 1/7 1/7 1/7 1/7 1/7]; % size 7

b3 = [1/12 1/12 1/12 1/12 1/12 1/12 1/12 1/12 1/12 1/12 1/12 1/12]; % size 12

y1 = filter(b1, a, x);
y2 = filter(b2, a, x);
y3 = filter(b3, a, x);

figure;
plot(t,x,t,y1);
title('Noised and Filtered signal (window 5)');
legend('Noised Signal', 'Filtered signal');

figure;
plot(t,x,t,y2);
title('Noised and Filtered signal (window 7)');
legend('Noised Signal', 'Filtered signal');

figure;
plot(t,x,t,y3);
title('Noised and Filtered signal (window 12)');
legend('Noised Signal', 'Filtered signal');
