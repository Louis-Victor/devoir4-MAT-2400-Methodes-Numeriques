f=@(x) 1/(1+x^2);
domain = -5:0.01:5;

plot(domain,f(domain))
hold on;

for n = [2, 4, 8]
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
hold off;
