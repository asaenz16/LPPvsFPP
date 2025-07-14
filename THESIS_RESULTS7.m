% Joint density of (t_min, t_max) for Brownian motion on [0,1]
% f(x1,x2) = (4/pi^2) * sum_{k1,k2,k3 >=1}
% [ (-1)^(k2+1) * k2^2 * (1-(-1)^k1)*(1-(-1)^k3) ]
% / [ k1^2 x1 + k2^2 (x2-x1) + k3^2 (1 - x2) ], x1 <= x2
% and the analogous swap when x1 > x2.
T = 1;
N = 1000; % truncate triple sum
A = 4/pi^2; % prefactor
% grid in [0,1]x[0,1]
x1_vals = linspace(0,T,100);
x2_vals = linspace(0,T,100);
[x1g, x2g] = meshgrid(x1_vals, x2_vals);
f_vals = zeros(size(x1g));
for i = 1:numel(x1g)
x1 = x1g(i);
x2 = x2g(i);
% if x1 > x2, swap roles of x1,x2 in the denominator
if x1 <= x2
denom_fun = @(n1,n2,n3) (n1^2*x1 + n2^2*(x2-x1) + n3^2*(1 - x2))^2;
else
denom_fun = @(n1,n2,n3) (n1^2*x2 + n2^2*(x1-x2) + n3^2*(1 - x1))^2;
end
s = 0;
for n1 = 1:N
if mod(n1,2)==0, continue; end % kills even k1
for n2 = 1:N % keep all k2
for n3 = 1:N
if mod(n3,2)==0, continue; end % kills even k3
denom = denom_fun(n1,n2,n3);
num = (-1)^(n2+1) * n2^2 * 4;
% the *4* is (1-(-1)^n1)*(1-(-1)^n3)=2*2 for odd n1,n3
s = s + num/denom;
end
end
end
f_vals(i) = A * s;
end
% plot
figure;
surf(x1g, x2g, f_vals, 'EdgeColor','none');
xlabel('t_{min} / T');
ylabel('t_{max} / T');
zlabel('f(t_{min},t_{max})');
title('Joint density of argmin and argmax (T=1)');
view(45,30);
colorbar;
axis([0 1 0 1 -1 5]);