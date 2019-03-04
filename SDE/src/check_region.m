function region = check_region(x,regions,dist)
    num_regions = size(regions,2);
    found_region = false;
    region = 0; %default is no region
    i = 1;
    while i <= num_regions && ~found_region
        if norm(x - regions(:,i)) < dist
            found_region = true;
            region = i;
        else
            i = i+1;
        end
    end
end