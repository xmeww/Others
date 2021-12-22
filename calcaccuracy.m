function out = calcaccuracy(X)

levels = [-10 -3.33 3.33 10];
nLevels = numel(levels);
out = zeros(1,nLevels);

for i = 1:nLevels
    out(i) = sum(X == levels(i))/numel(X);
end

end

