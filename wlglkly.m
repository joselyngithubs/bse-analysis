function negLogL = wlglkly(guess)

global val k n minimum max_minus_min

A = guess(1);
B = guess(2);
if A <= 0 || B <= 0
    negLogL = inf;
    return;
end

P = minimum + max_minus_min*(1-exp(-(val/A).^B));

negLogL = -log(P(:))'*k(:) - log(1-P(:)')*n(:);

end
