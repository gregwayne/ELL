function cst = mixing_loss(e,m,p, all_terms_meas)

    [calc] = mixing_all_terms(e,m,p);    
    nT          = sum(all_terms_meas);    
    cst         = sum((all_terms_meas-nT*calc).^2./abs(nT*calc));

end
