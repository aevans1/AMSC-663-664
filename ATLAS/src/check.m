a = rand(2,9)';
b = rand(2,9)';


% Check input dimensionality
[na,aDims] = size(a);
[nb,bDims] = size(b);
if (aDims ~= bDims)
    error('Input matrices must have the same dimensionality.');
end

% Create index matrices
[i,j] = meshgrid(1:na,1:nb);

% Compute array of inter-point differences
delta = a(i,:) - b(j,:);

% Compute distance by specified method
dmat = zeros(na,nb);

dmat = sqrt(sum(delta.^2,2))