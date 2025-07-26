function y = MSTHESIS_OPT_BARRIER(H,f,A,bd,x0)
P=H;
q=f;
MAXITERS =  100;
BETA = 0.65;
t = 4;
MU=5;
x=x0';
RESTOL=1e-18;
while(t < 1000)
    for iter = 1:MAXITERS
        j = bd-A*x;
        grad = t*(P*x+q) + A'*(1./j);
        hess = t*P + A'*diag(1./j.^2)*A;
        v = -hess\grad;
        s = 1; dj = -A*v;
        if(norm(v)<RESTOL)
            break;
        end
        while (min(j+s*dj) <= 0)
            s = BETA*s;
        end
        x = x+s*v;
        
    end
    t = MU*t;
%     disp(x)
end
% disp('Feasible Condition is Found')
y = x;
end