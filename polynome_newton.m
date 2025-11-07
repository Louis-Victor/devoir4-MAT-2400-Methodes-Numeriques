    function [c, a, D, fx] = polynome_newton(xi, fxi, x, resul)
%
%  Polynôme d'interpolation de Newton     
%  Programmeur: A. Lioret 
%  Référence: Analyse numérique pour ingénieurs, A. Fortin,
%             Presses internationales Polytechnique, 
%             Cinquième édition, 2015
%             Section 5.4 
%
%  function [c, a, D, fx] = polynome_newton(xi, fxi, x, resul)
%  Arguments
% 
%  Entree:
%  1) xi: Un vecteur ligne contenant les abscisses xi.
%     Exemple: [0 1 2 3].
%  2) fxi: Un vecteur ligne contenant les valeurs f(xi)
%     de la fonction à interpoler pour les xi correspondants. 
%     Exemple: [1 2 9 28].
%  3) x: Si vous voulez évaluer le polynôme de Newton en certains points,
%     vous devez fournir ces points (sinon, tapez [] pour avoir une matrice
%     vide). Exemple: [2.2 8.4].
%  4) resul: (facultatif) Si vous voulez les résultats dans un fichier, 
%     vous devez fournir le nom entre apostrophes (' ') de ce fichier.
%     Exemple: 'resul.dat' qui créera un fichier nommé resul.dat.
%  Sortie:
%  1) c: Vecteur contenant les coefficients du polynôme d'interpolation.
%     Exemple: c = [c3 c2 c1 c0] tel que
%                p3(x) = c3*x^3 + c2*x^2 + c1*x + c0.
%  2) a: Vecteur contenant les coefficients du polynôme  de Newton.
%     Exemple: a = [a0 a1 a2 a3] tel que
%                p3(x) =  a0
%                       + a1*(x-x0)
%                       + a2*(x-x0)*(x-x1)
%                       + a3*(x-x0)*(x-x1)*(x-x2)
%  3) D: Matrice contenant la table de différences divisées.
%  4) fx: Si on a fourni une vecteur x à évaluer, fx est le vecteur
%     contenant les évaluations du polynôme de Newton en ces 
%     points.
%
%  Exemples d'appel:
%  [c, a, D] = polynome_newton([0 1 2 3], [1 2 9 28], [])
%  [c, a, D] = polynome_newton([0 1 2 3], [1 2 9 28], [], 'resul.dat')
%  [c, a, D, fx] = polynome_newton([0 1 2 3], [1 2 9 28], [2.2 8.4])
%  [c, a, D, fx] = polynome_newton([0 1 2 3], [1 2 9 28], [2.2 8.4], 'resul.dat')

    n = length(xi);
    D = zeros(n,n);
    a = zeros(1,n);

    D(:,1) = fxi';
    a(1) = D(1,1);
    for j=2:n
       for k=2:j
          D(j,k) = (D(j,k-1)-D(j-1,k-1))/(xi(j)-xi(j-k+1));
       end
       a(j) = D(j,j); 
    end 

    c = D(n,n);
    for k=(n-1):-1:1
       c = conv(c,poly(xi(k)));
       h = length(c);
       c(h) = c(h) + D(k,k);
    end

    if ~(isempty(x))
       m = length(x);
       ev = zeros(n,m); 
       for l=1:m
          ev(1,l) = D(1,1);
          w = 1;       
          for k=2:n
             w = (x(l) - xi(k-1))*w;
             ev(k,l) = ev(k-1,l) + a(k)*w;
          end    
       end 
       fx = ev(n,:);
    end

    if (nargin > 3)
      if ~(isempty(x))
        polynome_newton_out(resul, xi, fxi, D, a, x, ev);
      else
        polynome_newton_out(resul, xi, fxi, D, a);
      end
    end
