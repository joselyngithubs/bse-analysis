function threshold = PitchCompareWeibull(Data)

DataDim = size(Data);
MinProb = .5;
MaxProb = .98;
Incorrect = Data(:,3)~=Data(:,4);
Correct = Data(:,3)==Data(:,4);
AbsCents = abs(log2(Data(:,2)./Data(:,1)))*1200;
WeibullData = [AbsCents Correct Incorrect];
[A, B] = FitWeibull(WeibullData,MinProb,MaxProb);

threshold = plotWeibullFunction(WeibullData,A,B);

end

function threshold = plotWeibullFunction(WeibullData,A,B)

[~,index] = sort(WeibullData(:,1));
sorteddata = WeibullData(index,:);

nBins = 7;
means = NaN(nBins,1);
propCorrect = NaN(nBins,1);

for i=1:nBins % if 7 bins, 10 trials per bin
    group = sorteddata((i-1)*10+1:i*10,:);
    means(i) = mean(group(:,1));
    propCorrect(i) = sum(group(:,2))/(sum(group(:,2))+sum(group(:,3)));
end

domain = 0:0.01:max(means);
f = Weibull(domain,A,B,0.5,0.98);

I = find(f<.8,1,'last'); % find index of 80% threshold
if isempty(I)
    threshold = NaN;
else
    threshold = domain(I);
end

end

function [A, B, A_firstguess, B_firstguess, NegLogLhood] = FitWeibull(Data,minprob,maxprob,firstGuessA)

% [A, B] = FitWeibull(data,minprob,maxprob)
%
% Obtain maximum likelihood fit of the Weibull for 2-alternative forced 
% choice data stored in array, data. The rows of data correspond to the set of 
% experimental conditions. Each condition is marked by 
% a stimulus value (e.g., a contrast level, a percent of stimulus dots, etc.).
% Each row of data comprises 3 columns. The first column gives 
% the stimulus value of the condition; the 2nd column gives the number 
% of correct responses observed for this condition, and the 3rd column gives 
% the number of incorrect responses observed for this condition. 

global val k n minimum max_minus_min

if maxprob == 1
    error('maxprob must be less than 1: suggested value: 0.98')
end
k = Data(:,2);
n = Data(:,3);
NumObservations = k+n;
val = Data(NumObservations > 0,1);
k = Data(NumObservations > 0,2);
n = Data(NumObservations > 0,3);
Trials = k+n;
NumConditions = length(val);
minimum = minprob;
max_minus_min = maxprob-minprob;

% We first obtain an approximation of A and B using a least-squares method.
p = k./Trials;
p_app1 = (maxprob - .001)*(p >= maxprob);
p_app0 = (minprob + .001)*(p <= minprob);
p_rest = p.*((p < maxprob) & (p > minprob));
p = p_app1 + p_app0 + p_rest;
y = log(-log(1-(p-minimum)/max_minus_min));
X = [ones(length(val),1), log(val)];
% b = regress(y,X); % This line requires statistics toolbox.  If you don't have
	% it, use the following:
b = pinv(X)*y;

B_firstguess = b(2);
if nargin < 4
    A_firstguess = exp(-(b(1)/b(2)));
else
    A_firstguess = firstGuessA;
end
first_guess = [A_firstguess, B_firstguess];

OptStruct = optimset('MaxFunEvals',100000,'MaxIter',10000,'TolFun',0.0001,'TolX',0.0001);
[A_B,NegLogLhood] = fminsearch('wlglkly', first_guess,OptStruct);
A = A_B(1);
B = A_B(2);

end

function y = Weibull(x,A,B,minProb,maxProb)

y = minProb + (maxProb-minProb)*(1-exp(-(x/A).^B));

end