function result = logit(p)
% logit function: natural log(p/(1-p))

result = NaN(size(p));
for i=1:length(p)
    % corrections
    if p(i)==0
        p(i) = 0.01;
    end
    if p(i)==1
        p(i) = 0.99;
    end
    
    result(i) = log(p(i)/(1-p(i)));
end 

end