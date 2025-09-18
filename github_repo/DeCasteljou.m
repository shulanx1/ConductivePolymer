
function [ct] = DeCasteljou(P, t)

n = size(P,2);
for i=1:n
   for j=1:n-i
       P(:,j) = (1-t)*P(:,j) + t*P(:,j+1);
   end
end
ct = P(:,1);

end

