load('current_atlas.mat');
luke_b = new_S.b;
luke_c = new_S.c;
luke_sig = new_S.sigma;
luke_T = new_S.T;
luke_mu = new_S.mu;

load('LearnedSimulator_CM2D.mat');


for j = 1:218
   
    if norm(b(:,j) - luke_b(:,j)) > 0.000000001
        j
    end
    
end

for i = 1:218
    for j = 1:218
        if norm(c(:,i,j) - luke_c(:,i,j)) > 0.000000001
            fprintf("(%d,%d) \n",i,j)
        end
    end
end


for i = 1:218
    for j = 1:218
        if norm(mu(:,i,j) - luke_mu(:,i,j)) > 0.000000001
            fprintf("(%d,%d) \n",i,j)
        end
    end
end



for i = 1:218
        if norm(sig(:,:,i) - luke_sig(:,:,i)) > 0.000000001
            fprintf("(%d,%d) \n",i)
        end
end

for i = 1:218
    for j = i + 1:218
        if norm(T(:,:,i,j) - luke_T(:,:,i,j)) > 0.000001
            fprintf("(%d,%d) \n",i,j)
        end
    end
end



