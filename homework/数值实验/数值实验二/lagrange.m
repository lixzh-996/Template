%----------------------plot origin function---------------------%
a=-1:0.001:1;
b=1./(1+25.*a.^2);
plot(a,b,'k','LineWidth',1.5);
hold on
%---------------Lagrange interpolating polynomial---------------%
n=[1,2,3,4,5,6,7]; %nodes%
for i = n
    x=-1:2/i:1;
    y=1./(1+25.*x.^2);
    x_2=-1:0.05:1;
    lagrange(x,y,x_2);
end
legend('origin','n=2','n=3','n=4','n=5','n=6','n=7')
%----------------------constructe polynomial--------------------%
function L = lagrange(x,y,x_2)
a=x_2;
L=zeros(1,length(a));

for i = 1:length(a)
    l = ones(1,length(x));
    for k = 1:length(x)
        for j = 1:length(x)
            if j ~= k
                l(k)=l(k)*(a(i)-x(j))/(x(k)-x(j));
            end
        end
        L(i)=L(i)+l(k)*y(k);
    end
end

plot(a,L,'--');
end






