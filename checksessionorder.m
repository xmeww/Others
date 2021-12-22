function out  = checksessionorder(sessionOrder)
% Checks if all the required sessions in the required numbers are present
out = false(size(sessionOrder,1),1);

for i = 1:size(sessionOrder,1)
    template = [1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,6,7,8];
    
    out(i) = all(sort(sessionOrder(i,:)) == template);
end


end

