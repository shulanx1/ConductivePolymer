

n = 400;

P = [1 2 3 4 5 6;
    8 1 3 5 1 9];

X = [];
t = 0;
tIncr = 1/400;
for i=1:n
    
    X(:,i) = DeCasteljou(P,t);
    t = i*tIncr;
    
end

plot(P(1,:),P(2,:),'-or',X(1,:),X(2,:),'-b');