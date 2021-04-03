function dprime = analyzeDprime(data)

types = data(:,1);
resps = data(:,2);

nSignalTrials = sum(types==2); % major
nNoiseTrials = sum(types==1); % minor

nHits = sum(types==2 & resps==2);
nMisses = sum(types == 2 & resps == 1);
nCRs = sum(types==1 & resps==1);
nFAs = sum(types==1 & resps==2);

pHit = nHits/nSignalTrials;
pMiss = nMisses/nSignalTrials;
pCR = nCRs/nNoiseTrials;
pFA = nFAs/nNoiseTrials;

if pHit == 1
    pHit = (nSignalTrials - 0.5) / nSignalTrials;
end
if pFA == 1
    pFA = (nNoiseTrials - 0.5) / nNoiseTrials;
end
if pHit == 0
    pHit = 0.5 / nSignalTrials;
end
if pFA == 0
    pFA = 0.5 / nNoiseTrials;
end

dprime = norminv(pHit)-norminv(pFA);

end