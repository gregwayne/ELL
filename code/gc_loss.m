function cst = gc_loss(e,m,p)

    n       = 1 - e - m - p;
    
    nE      = 101;
    nM      = 4;
    nP      = 7;
    nEM     = 20;
    nEP     = 12;
    nPM     = 0;
    nEMP    = 0;
    nN      = 1;
    
    nT      = nE+nM+nP+nEM+nEP+nPM+nEMP+nN;
   
    unival  = @(x,y) (x^3 + 3*x^2*y + 3*x*y^2);
    bival   = @(x,y,z) (3*x^2*y + 3*x*y^2 + 6*x*y*z);
    trival  = @(x,y,z) x*y*z;
    
    cst     =   (nE  - nT*unival(e,n))^2 ...
              + (nM  - nT*unival(m,n))^2 ...
              + (nP  - nT*unival(p,n))^2 ...
              + (nEM - nT*bival(e,m,n))^2 ...
              + (nEP - nT*bival(e,p,n))^2 ...
              + (nPM - nT*bival(p,m,n))^2 ...
              + (nEMP - nT*trival(e,m,p))^2 ...
              + (nN  - nT*trival(n,n,n))^2;

end