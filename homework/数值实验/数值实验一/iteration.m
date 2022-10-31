%------------------------------Iteration----------------------------%
%-------------------------------------------------------------------%
%-------------------------parameter settings------------------------%
h = input('Enter the step h:');
N = 1/h;
e = input('Enter the error e:');

%------------------------inital Matrix L_2--------------------------%
C = zeros(N-1);
for i = 1:N-1
    for j = 1:N-1
        if abs(i-j)==1
            C(i,j)=1;
        end
    end
end
I = eye(N-1);

L = zeros((N-1)*(N-1));
row = (N-1)*ones(1,N-1);
L_1 = mat2cell(L,row,row);

for i = 1:N-1
    L_1{i,i} = 4*I-C;
end

for i = 1:N-1
    for j = 1:N-1
        if abs(i-j)==1
            L_1{i,j}=-I;
        end
    end
end

L_2 = cell2mat(L_1);

%--------------------------initial x,y,f----------------------------%
x = zeros(1,N-1);
y = zeros(1,N-1);

for i = 1:N-1
    x(i)=i*h;
    y(i)=i*h;
end

f = zeros((N-1)^2,1);
for j = 1:N-1
    for i = 1:N-1
        f((j-1)*(N-1)+i) = h^2*2*pi^2*sin(pi*x(i))*sin(pi*y(j));
    end
end

%--------------------------initial vetor----------------------------%
u_0 = zeros((N-1)^2,1);  %initial vector%
u_e = zeros((N-1)^2,1);  %exact solution vector%
for j = 1:(N-1)
    for i = 1:(N-1)
        u_e((j-1)*(N-1)+i) = sin(pi*x(i))*sin(pi*y(j));
    end
end

%-------------------------matrix D,L,U------------------------------%
D = zeros((N-1)^2,(N-1)^2);
for i = 1: (N-1)^2
    D(i,i) = L_2(i,i);
end

L = zeros((N-1)^2,(N-1)^2);
for i = 1:(N-1)^2
    for j = 1:(N-1)^2
        if j<i
            L(i,j) = -L_2(i,j);
        end
    end
end

U = zeros((N-1)^2,(N-1)^2);
for i = 1:(N-1)^2
    for j = 1:(N-1)^2
        if j>i
            U(i,j) = -L_2(i,j);
        end
    end
end

%------------------------Jacobi iteration---------------------------%
J = pinv(D)*(L+U);
f_j = pinv(D)*f;
u_1 = u_0;
u_2 = f_j;
n = 1;   %  iteration number  %
while get_norm(u_2,u_1)>=e
     u_1 = J*u_2 + f_j;
     t = u_2;
     u_2 = u_1;
     u_1 = t;
     n = n+1;
end

rhoJ = max(abs(eig(J)));
n_J = n;
nm = get_norm(u_2,u_e);   %||u^(k+1)-u*||%
fprintf('the Jacobi iteration number is %d',n_J);
fprintf('the norm of the error is %f',nm);
fprintf('the radius of the convergence is %f',rhoJ);

%----------------------Gauss-Seidel iteration-----------------------%
G = pinv(D-L)*U;
f_g = pinv(D-L)*f;
u_1 = u_0;
u_2 = f_g;
n = 1;   %  iteration number  %
while get_norm(u_2,u_1)>=e
     u_1 = G*u_2 + f_g;
     t = u_2;
     u_2 = u_1;
     u_1 = t;
     n = n+1;
end

rhoG = max(abs(eig(G)));
n_G = n;
nm = get_norm(u_2,u_e);   %||u^(k+1)-u*||%
fprintf('the Gauss-Seidel iteration number is %d',n_G);
fprintf('the norm of the error is %f',nm);
fprintf('the radius of the convergence is %f',rhoG);

%-------------------------------SOR---------------------------------%
omega = input('Enter the weight \omega:');
L_omega = pinv(D-omega*L)*((1-omega)*D+omega*U);
f_sor = omega*pinv((D-omega*L))*f;
u_1 = u_0;
u_2 = f_sor;
n = 1;   %  iteration number  %
while get_norm(u_2,u_1)>=e
     u_1 = L_omega*u_2 + f_sor;
     t = u_2;
     u_2 = u_1;
     u_1 = t;
     n = n+1;
end

rhoS = max(abs(eig(L_omega)));
n_SOR = n;
nm = get_norm(u_2,u_e);   %||u^(k+1)-u*||%
fprintf('the SOR iteration number is %d',n_SOR);
fprintf('the norm of the error is %f',nm);
fprintf('the radius of the convergence is %f',rhoS);
%---------------------------function--------------------------------%
function norm = get_norm(a,b)
    norm = max(abs(a-b));
end










    