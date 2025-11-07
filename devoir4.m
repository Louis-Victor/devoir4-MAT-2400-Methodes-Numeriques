f=@(x) 1/(1+x^2);
domain = -5:0.01:5;

fedex=arrayfun(@(x) f(x),domain);


plot(domain,fedex)
hold on;

N = [2,4,8];
nomLeg = arrayfun(@(k) num2str(k),N);

for n=N
    n
    h=10/n;

    xi = [];
    fxi = [];
    for i=0:n
        xi(i+1) = -5+i*h;
        fxi(i+1) = f(xi(i+1));
    end

    [a,b,c,utile] = polynome_newton(xi,fxi,domain);

    plot(domain,utile);
    
end
legend("f(x)",["P_",nomLeg(1)],["P_",nomLeg(2)],["P_",nomLeg(3)])
hold off;
