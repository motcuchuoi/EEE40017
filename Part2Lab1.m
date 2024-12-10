%step = zeros([1,49]);
%step(10:49) = 1;

%impulse = zeros([1, 49]);
%impulse(10) = 1;

% causal system y(n) = 0 with n < 0 
% y(n) = y(n-1) + 1/k(x(n)-x(n-k));
% Unit sample u(n) = 1 (aka x(n) for n >= 0)

% plot the output
%Y = filter([1,1,1,1,1],5,impulse);

%stem(1:49,Y,"LineWidth",2);

% values + range
n = (-10:39);
k = 5;
for m = (0:5): % voi i chay tu 0 den 5 va i thuoc tap n(-10:39)
    y(m,k)
 
a = y(n,k);
if n > 0:
    a = 1/k(x(n)+x(n-1)+)
