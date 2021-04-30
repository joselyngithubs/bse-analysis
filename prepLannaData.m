function compressedData = prepLannaData

% read in one reasonable subject's data to be used in MCMC
dataLoc = 'C:\Users\Joselyn\Documents\GitHub\BSE-data\mat\';
files = dir([dataLoc '*.mat']);
files = {files.name};

nSubj = length(files);

lanna_all = [];

f=2; % subject 2
load([dataLoc files{f}]);
lanna_all = [lanna{3,2},lanna{2,2}==1,lanna{2,2}==2];

compressedData = NaN(6,3);
for k=1:6
    whichOnes = lanna_all(:,1)==k;
    compressedData(k,1) = k; % col 1: stimulus type (1 thru 6)
    compressedData(k,2) = sum(lanna_all(whichOnes,3)); % col 2: num trials responded "2"
    compressedData(k,3) = sum(lanna_all(whichOnes,2)); % col 3: num trials responded "1"
end

end