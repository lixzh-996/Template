%----------------------------initial parameter----------------%
L = 2;
W = 1;
Nx = input('Enter the Nx:'); 
Ny = input('Enter the Ny:'); 

h_x = L/Nx;
h_y = W/Ny;

%----------------------------initial Matrix-------------------%
%initial C%
C = zeros(Nx-1); 
for i = 1:Nx-1 
    C(i,i)=-2*(1/(h_x)^2+1/(h_y)^2);
    for j = 1:Nx-1
        if abs(i-j)==1
            C(i,j)=1/(h_x)^2;
        end
    end
end

%initial A%
I = eye(Nx-1);

B = zeros((Nx-1)*(Ny-1));
row = (Nx-1)*ones(1,Ny-1);
B_1 = mat2cell(B,row,row);

for i = 1:Ny-1
    B_1{i,i} = C;
end

for i = 1:Ny-1
    for j = 1:Ny-1
        if abs(i-j)==1
            B_1{i,j}=1/(h_y)^2*I;
        end
    end
end

A = cell2mat(B_1);

%-----------------initial x,y,b-----------------%
x = zeros(1,Nx-1);
y = zeros(1,Ny-1);
phi_0=0.1;

for i = 1:Nx-1
    x(i)=i*h_x;
end

for j = 1:Ny-1
    y(j)=j*h_y;
end

%initial b%
b = zeros((Nx-1)*(Ny-1),1);
for i = 1:Nx-1
    for j = 1:Ny-1
        if i == Nx-1
            b((j-1)*(Nx-1)+Nx-1) = b((j-1)*(Nx-1)+Nx-1)+1/(h_x)^2*phi_0;
        end
        if j == Ny-1
            b((j-1)*(Nx-1)+i) = b((j-1)*(Nx-1)+i)+i/(h_y^2*Nx)*phi_0;
        end
        if j == 1
            b((j-1)*(Nx-1)+i) = b((j-1)*(Nx-1)+i)+i/(h_y^2*Nx)*phi_0;
        end
    end
end
b=-b;
%----------------solve the equation------------%
%--------------------\rho=0--------------------%
phi = A\b;

%---------------------\rho---------------------%
%--------------set e_0 =1,rho_0=10--------------%
rho = zeros((Nx-1)*(Ny-1),1);
for i = 1:Nx-1
    for j = 1:Ny-1
        rho((j-1)*(Nx-1)+i) = 10*exp(-20*((i*h_x-L/2)^2+(j*h_y-W/2)^2)/W^2);
    end
end

phi_1 = A\(b-rho);

%---------------------plot---------------------%
map1 = reshape(phi',[Nx-1,Ny-1])';
map2 = reshape(phi_1',[Nx-1,Ny-1])';

[X, Y] = meshgrid(x, y);
surf(X,Y,map1);
surf(X,Y,map2);

