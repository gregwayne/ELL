function cst = mixing_loss_4types(e,m,l,p, all_terms_meas)

    [calc] = mixing_all_terms_4types(e,m,l,p);    
    nT          = sum(all_terms_meas);    
    cst         = sum((all_terms_meas-nT*calc).^2./abs(nT*calc));

end
