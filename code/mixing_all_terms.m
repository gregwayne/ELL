function [all_terms_calc] = mixing_all_terms(e,m,p)

    n       = 1 - e - m - p;

    unival  = @(x,y) (x^3 + 3*x^2*y + 3*x*y^2);
    bival   = @(x,y,z) (3*x^2*y + 3*x*y^2 + 6*x*y*z);
    trival  = @(x,y,z) 6*x*y*z;

    all_terms_calc      = [unival(e,n), unival(m,n), unival(p,n) ...
                        bival(e,m,n), bival(e,p,n), bival(p,m,n) ...
                        trival(e,m,p), n^3];
    
end
