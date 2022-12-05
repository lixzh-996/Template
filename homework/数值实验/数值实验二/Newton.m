%----------------------Newton iteration----------------%
x_0=zeros(3,1);
f=-Vect(x_0);
F=mat(x_0);
dx=F\f;
x_1=x_0+dx;
n=0;
while norm(x_1-x_0)>=10^(-5)
    Dx=mat(x_1)\(-Vect(x_1));
    t=x_1;
    x_1=Dx+t;
    x_0=t;
    n=n+1;
end

function F = mat(x)

F=ones(3,3);
F(1,1)=3;
F(1,2)=x(3)*sin(x(2)*x(3));
F(1,3)=x(2)*sin(x(2)*x(3));
F(2,1)=2*x(1);
F(2,2)=-81;
F(2,3)=-cos(x(3));
F(3,1)=-x(2)*exp(-x(1)*x(2));
F(3,2)=-x(1)*exp(-x(1)*x(2));
F(3,3)=20;

end

function f = Vect(x)

f=ones(3,1);
f(1)=3*x(1)-cos(x(2)*x(3))-0.5;
f(2)=x(1)^2-81*(x(2)+0.1)+sin(x(3))+1.06;
f(3)=exp(-x(1)*x(2))+20*x(3)+(10*pi-3)/3;

end



















