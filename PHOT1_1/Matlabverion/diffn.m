function f = diffnn(x)

n = length(x);
f(1:n) = 0;

f(2:n-1) = 0.5*(x(3:n)-x(1:n-2));
f(1) = x(2)-x(1);
f(n) = x(n)-x(n-1);
% f(1:n-1) = diff(x);
% f(n) = x(n)-x(n-1);
f = f';
return;

