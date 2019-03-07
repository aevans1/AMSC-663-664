
    %Build array of wells assigned to charts
    wells = zeros(1,size(net,2));
    for n = 1:size(net,2)
        wells(n) = check_region(net(:,n),regions,dist);
    end