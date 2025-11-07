    function [S, fpp, Sx] = spline_cubique(xi, fxi, cdt, x, resul)
%
%  Interpolation par les splines cubiques    
%  Programmeur: A. Fortin 
%  Référence: Analyse numérique pour ingénieurs, A. Fortin,
%             Presses internationales Polytechnique, 
%             Cinquième édition, 2015
%             Section 5.6 
%
%  function [S, fpp, Sx] = spline_cubique(xi, fxi, cdt, x, resul)
%  Arguments
%
%  Entrée:
%  1) xi: Un vecteur ligne contenant les abscisses xi.
%     Exemple: [1 2 4 5].
%  2) fxi: Un vecteur ligne contenant les valeurs f(xi)
%     de la fonction à interpoler pour les xi correspondants. 
%     Exemple: [1 4 9 11].
%  3) cdt: Le type de spline cubique, soit:
%              'nat'       pour la spline naturelle (f0"=fn"=0)
%              'der_prem'  pour la spline telle que f0'=a et fn'=b
%              'der_sec'   pour la spline telle que f0"=a et fn"=b
%              'courb_cte' pour la spline telle que f0"=f1" et fn-1"=fn" 
%                          (courbure constante)
%              'not'       pour la condition not-a-knot
%              'per'       pour une spline périodique
%  4) x: Si vous voulez évaluer la spline cubique en certains points,
%     vous devez fournir ces points (sinon, tapez [] pour avoir une
%     matrice vide). Exemple: [2.2 8.4].
%  5) resul: (facultatif) Si vous voulez les résultats dans un fichier, 
%     vous devez fournir le nom entre apostrophes (' ') de ce fichier.
%     Exemple: 'resul.dat' qui créera un fichier nommé resul.dat.
%  Sortie:
%  1) S: Matrice contenant les coefficients des 
%     polynômes constituant la spline.
%     La i-ième ligne de S est S(i,:) = [fi  fi' fi'' fi''']
%     et telle que 
%     pi(x) = fi + fi'(x-xi) + (fi''/2)(x-xi)^2 + (fi'''/6)(x-xi)^3     
%  2) fpp:  Vecteur contenant les dérivées secondes de 
%     la spline cubique.
%  3) Si on a fourni un vecteur x où évaluer, Sx est le vecteur
%     contenant les évaluations de la spline cubique en ces 
%     points.
%
%  Exemples d'appel:
%  [S, fpp]     = spline_cubique([1 2 4 5], [1 4 9 11], 'nat', [])
%  [S, fpp]     = spline_cubique([1 2 4 5], [1 4 9 11], 'nat', [], 'resul.dat')
%  [S, fpp, Sx] = spline_cubique([1 2 4 5], [1 4 9 11], 'nat', [2.2 3.18])
%  [S, fpp, Sx] = spline_cubique([1 2 4 5], [1 4 9 11], 'nat', [2.2 3.18], 'resul.dat')
tic
    n = length(xi);
    [xi,ind] = sort(xi);
    fxi = fxi(ind);

    if strcmp(cdt, 'der_prem') == 1
      Mx0 = 'Entrez f0'' et appuyez sur la touche ENTER. f0'' = ';
      fpx0 = input(Mx0);
      Mxn = ['Entrez f',int2str(n),''' et appuyez sur la touche ENTER. f', int2str(n), ''' = '];
      fpxn = input(Mxn);
    end
    if strcmp(cdt, 'der_sec') == 1
      Mx0 = 'Entrez f0" et appuyez sur la touche ENTER. f0" = ';
      fppx0 = input(Mx0); 
      Mxn = ['Entrez f',int2str(n),'" et appuyez sur la touche ENTER. f', int2str(n), '" = '];
      fppxn = input(Mxn);
    end
%   Calcul des h_i
    h = diff(xi);
    
%   Première différence divisée de f
    D = diff(fxi)./h;   
%   Calcul de la tridiagonale
    diviseur = (h(1:n-2) + h(2:n-1));
%   Diagonle inférieure
    diag_inf = h(1:n-2)./diviseur;
%   Diagonale supérieure
    diag_sup = h(2:n-1)./diviseur;
%   Terme de droite (différences divisées d ordre 2)
    V = 6*diff(D)./diviseur;

    A = sparse(n,n);
  %  A = zeros(n,n);
    b = zeros(n,1);
    for i=2:n-1
	   A(i,i) = 2.0;
	   A(i,i-1)= diag_inf(i-1);
	   A(i,i+1) = diag_sup(i-1);
	   b(i) = V(i-1); 
    end
   
    if strcmp(cdt, 'nat') == 1
         A(1,1) =1;
	     A(n,n) = 1;
         b(1) = 0;
         b(n) = 0;
    elseif strcmp(cdt, 'der_prem') == 1
	    A(1,1) = 2;
	    A(1,2) = 1;
	    b(1) = (6/h(1))*( D(1) - fpx0);
	    A(n,n) = 2;
	    A(n,n-1) = 1;
	    b(n) = (6/h(n-1))* (fpxn - D(n-1));
    elseif  strcmp(cdt, 'courb_cte') == 1
	    A(1,1) = 1;
        A(1,2) = -1;
	    A(n,n) = 1;
	    A(n,n-1) = -1;
	    b(1) = 0;
        b(n) = 0;
    elseif  strcmp(cdt, 'der_sec') == 1
        A(1,1) =1;
	    A(n,n) = 1;
	    b(1) = fppx0;
	    b(n) = fppxn;
    elseif   strcmp(cdt, 'not') == 1
	    A(1,1) = h(2);
        A(1,2) = -diviseur(1);
        A(1,3) = h(1);
        A(n,n-2)= h(n-1);
        A(n,n-1) = -diviseur(n-2);
        A(n,n) = h(n-2);
        b(1) = 0;
	    b(n) = 0;
    elseif   strcmp(cdt, 'per') == 1
	    A(1,1) = 1;
        A(1,n) = -1;
	    A(n,1) = h(1)/3;
        A(n,2) = h(1)/6;
	    A(n,n-1) = h(n-1)/6;
        A(n,n) = h(n-1)/3;
        b(1) = 0;
	    b(n) = D(1) - D(n-1);
    else
        error('Mauvais nom pour les conditions aux limites');
    end

%   Résolution du système
%format rat
     % A
	  fpp = A\b;
      S = zeros(1);
%
      for k=1:n-1
          S(k,1) = fxi(k);
          S(k,2) = D(k) - h(k)*fpp(k)/3 - h(k)*fpp(k+1)/6;
          S(k,3) = fpp(k);
          S(k,4) = ((fpp(k+1)-fpp(k))/h(k));
      end
%    
%    Évaluation de la spline pour le vecteur x
    for k=1:n-1
      I = find((x<=xi(k+1)) & (x>=xi(k)));
      w = x(I) - xi(k);
      Sx(I) = S(k,1) + w.*(S(k,2) + w.* (S(k,3)/2. + w.* S(k,4)/6) ); % Horner
    %  Sx(I) = S(k,1) + w.*S(k,2) + w.^2 *S(k,3)/2. + w.^3* S(k,4)/6;
    end
    

    if (nargin > 4)
      if strcmp(cdt, 'dp') == 1
        spline_cubique_out(resul, S, fpp, Sx, xi, fxi, cdt, x, fpx0, fpxn);
      else
        spline_cubique_out(resul, S, fpp, Sx, xi, fxi, cdt, x);
      end
    end
