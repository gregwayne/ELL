function [all_terms_calc] = mixing_all_terms_4types(e,m,l,p)

    n       = 1 - e - m - l - p;

unival=@(x,y)         (x^3 + 3*x^2*y + 3*x*y^2);
bival = @(x,y,z)      (6*z*x*y + 3*x*y^2 + 3*x^2*y);
trival = @(x,y,z)     (6*x*y*z);

all_terms_calc = [unival(e,n), unival(m,n), unival(l,n), unival(p,n) ...
                  bival(e,m,n), bival(e,l,n), bival(e,p,n), ...
                  bival(m,l,n), bival(m,p,n), bival(l,p,n), ...
                  trival(e,m,l), trival(e,m,p), trival(e,l,p), trival(m,l,p), ...
                  n^3];
    
end