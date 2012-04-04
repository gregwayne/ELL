function [all_terms_calc,all_terms_meas] = gc_all_terms(e,m,p)

    nE      = 101;
    nM      = 4;
    nP      = 7;
    nEM     = 20;
    nEP     = 12;
    nPM     = 0;
    nEMP    = 0;
    nN      = 12;
    
    all_terms_meas      = [nE,nM,nP,nEM,nEP,nPM,nEMP,nN];

    n       = 1 - e - m - p;

    unival  = @(x,y) (x^3 + 3*x^2*y + 3*x*y^2);
    bival   = @(x,y,z) (3*x^2*y + 3*x*y^2 + 6*x*y*z);
    trival  = @(x,y,z) x*y*z;

    all_terms_calc      = [unival(e,n), unival(m,n), unival(p,n) ...
                        bival(e,m,n), bival(e,p,n), bival(p,m,n) ...
                        trival(e,m,p), trival(n,n,n)];
    
end