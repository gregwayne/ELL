n_proto = 5;
n_pca   = 10;

[prototypes,bfns] =compute_granule_prototypes(n_proto,n_pca);
figure(2);

% ann says that the units start at -0.025 seconds before EOCD and
% end at 0.2 seconds afterwards
for i=1:n_proto
    subplot(n_proto/2,2,i);
    plot(linspace(-.025,.2,4500),bfns(prototypes(i),:));    
    axis([-0.025 0.2 -10 20]);
end

    xlabel('seconds');
    ylabel('mV');