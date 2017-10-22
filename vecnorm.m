function Y = vecnorm(X)
[~,n] = size(X);
Y = zeros(n,1);
for i=1:n
    Y(i) = norm(X(:,i));
end