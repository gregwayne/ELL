function cst = gc_loss(e,m,p)

    [calc,meas] = gc_all_terms(e,m,p);    
    nT          = sum(meas);    
    cst         = sum((meas-nT*calc).^2./abs(nT*calc));

end