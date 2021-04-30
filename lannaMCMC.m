function [samples, LogL, Sigmas, sigmaScalar] = lannaMCMC(Data, NumSamples, Prevs, Sigmas, sigmaScalar, fixSampling)

global NumParams D SigmaScalar

UpdateBlockSz = 5000;
NumParams = 2;  % mean and std

if nargin < 6 || isempty(fixSampling)
    fixSampling = false;
end

if nargin < 5 || isempty(sigmaScalar)
    SigmaScalar = 0.01;
end

if nargin < 4 || isempty(Sigmas)
    Sigmas = [0.01 0.01];
end

if nargin < 3 || isempty(Prevs)
    Prevs = [3.5 2];
end

if nargin < 2 || isempty(NumSamples)
    NumSamples = 100000;
end

Mins = [-10 0];
Maxes = [10 10];

D = Data;

samples = NaN(NumSamples,NumParams);
samples(1,:) = Prevs;
LogL(1,:) = GetLogLikelihood(Prevs);

UpdateBlockCtr = 0;
LRats = NaN(UpdateBlockSz,1);
for k = 2:NumSamples
    UpdateBlockCtr = UpdateBlockCtr+1;
    Candidate = GetCandidate(samples(k-1,:),Mins,Maxes,Sigmas);
    LogLhoodCand = GetLogLikelihood(Candidate);
    LRats(UpdateBlockCtr) = exp(LogLhoodCand - LogL(k-1));
    if rand < LRats(UpdateBlockCtr)   % keep the candidate sample
        samples(k,:) = Candidate;
        % LogLhoodCand
        LogL(k) = LogLhoodCand;
    else
        samples(k,:) = samples(k-1,:);
        LogL(k) = LogL(k-1);
    end
    if UpdateBlockCtr == UpdateBlockSz
        fprintf('iter = %d, median likelihood ratio = %0.2f\n',k,median(LRats));
        if ~fixSampling
            Sigmas = UpdateSigmas(Sigmas,samples(k-UpdateBlockSz+1:k,:),LRats);
        end
        UpdateBlockCtr = 0;
    end
end
sigmaScalar = SigmaScalar;
end

function LogL = GetLogLikelihood(Params)

global D

% D(:,1) = stimulus type (1 thru 6). 1-3=="1"; 4-6=="2"
% D(:,2) = num trials responded "2"
% D(:,3) = num trials responded "1"

% P(respond "2")
% = P(mu - sigma * X > 0)
% = P(mu > sigma * X)
% = P(X < mu/sigma)
% = normcdf(mu/sigma)

% P(respond "1")
% = 1 - normcdf(mu/sigma)

% Likelihood function
% = [(normcdf(mu/sigma) for all trials resp "2") * (1-normcdf(mu/sigma) for all trials resp "1")]

mu = Params(1);
sigma = Params(2);

prob = normcdf((D(:,2)-mu)/sigma); 
prob(prob > 0.99) = 0.99;
prob(prob < 0.01) = 0.01;
prob = prob(:);

LogL = sum(D(:,2).*log(prob) + D(:,3).*log(1-prob));

end

function    NextCandidate = GetCandidate(tmpC,Mins,Maxes,sigmas)

tmp = tmpC + randn(size(sigmas)).*sigmas;
cnt = 1;
while cnt < 100 && any(tmp < Mins) || any(tmp > Maxes)
    tmp = tmpC + randn(size(sigmas)).*sigmas;
    cnt=cnt+1;
end
if cnt >= 100
    error('cnt == 100');
end
NextCandidate=tmp;
end


function   Sigmas = UpdateSigmas(Sigmas,crntSamples,LRats)

global SigmaScalar

Med = median(LRats);
if Med < 0.35    % decrease scalar for sigma
    SigmaScalar = 0.9*SigmaScalar;
else
    SigmaScalar = SigmaScalar/0.9;
end
Sigmas = 0.5*Sigmas + 0.5*SigmaScalar*std(crntSamples);
end
