function [newX] = scaleFunc(X)

%scaleFact = 0.6*max(max(X(:)),-min(X(:)));  % another way to do it
%scaleFact = max(abs(X(:)));
scaleFact = prctile(abs(X(:)),95); %originally implemented way

newX = X ./ scaleFact;


end