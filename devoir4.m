% la fonction
f=@(x) 1/(1+x^2);
% le domainne de la fonction
domain = -5:0.01:5;

% l'évaluation de f sur le domaine
fedex=arrayfun(@(x) f(x),domain);


plot(domain,fedex)
hold on;

% le nombre de points pris pour l'interpolation
N = [2,4,8];
% nom des polynomes P_n pour la légende
%nomLeg = arrayfun(@(k) num2str(k),N);

% génération des interpolations
for n=N
    % équart entre les points
    h=10/n;

    % abscisse des points d'interpolation
    xi = [];
    % ordonée des points d'interpolation
    fxi = [];
    for i=0:n
        xi(i+1) = -5+i*h;
        fxi(i+1) = f(xi(i+1));
    end

    [a,b,c,utile] = polynome_newton(xi,fxi,domain);

    plot(domain,utile);
    
end
nomLeg = compose("P_%d",N);
legend(["f(x)",nomLeg],"Location","best")
%legend("f(x)",["P_",nomLeg(1)],["P_",nomLeg(2)],["P_",nomLeg(3)])
hold off;
