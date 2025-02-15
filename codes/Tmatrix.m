function T = Tmatrix(X)
% Gets the orientation matrix (T) from data in X.
% Parameters:
%   X : nested matrix of input data
% Returns:
%   T : orientation matrix as a nested matrix

T = zeros(3); % Initialize the orientation matrix

for i = 1:size(X, 1)
    row = X(i,:);
    for k = 1:3
        for l = 1:3
            T(k,l) = T(k,l) + row(k) * row(l);
        end
    end
end
end
