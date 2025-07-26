function y = MSTHESIS_OPT_IP_FULL(H,f,A,bd,x0)
P=H;
q=f;
[m,n]=size(A);

MAXITERS = 2000;
TOL = 1e-6;
RESTOL = 1e-8;
MU = 10;
ALPHA = 0.01;
BETA = 0.5;
x = x0'; s = bd-A*x; z = 1./s;
for iters = 1:MAXITERS
gap = s'*z; res = P*x + q + A'*z ;
if ((gap < TOL) && (norm(res) < RESTOL)), break; end
tinv = gap/(m*MU);
sol = -[ P A'; A diag(-s./z) ] \[ P*x+q+A'*z; -s + tinv*(1./z) ];
dx = sol(1:n); dz = sol(n+[1:m]); ds = -A*dx;
r = [P*x+q+A'*z; z.*s-tinv];
step = min(1.0, 0.99/max(-dz./z));
while (min(s+step*ds) <= 0)
    step = BETA*step; 
end

newz = z+step*dz; 
newx = x+step*dx; 
news = s+step*ds;
newr = [P*newx+q+A'*newz; newz.*news-tinv];

while (norm(newr) > (1-ALPHA*step)*norm(r))
step = BETA*step;
newz = z+step*dz; 
newx = x+step*dx; 
news = s+step*ds;
newr = [P*newx+q+A'*newz; newz.*news-tinv];
end
x = x+step*dx; z = z +step*dz; s = bd-A*x;
end

y=x;
end