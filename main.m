function main
% main analysis for BSE (conducted on mat files) which calls individual
% function files to analyze each task separately

dataLoc = 'C:\Users\Joselyn\Documents\GitHub\BSE-data\mat\';
files = dir([dataLoc '*.mat']);
files = {files.name};

nSubj = length(files);

% preallocate

pCorrAnchor = NaN(nSubj,1);
pCorrTest = NaN(nSubj,1);
pCorrAll = NaN(nSubj,1);

pCorrShape = NaN(nSubj,1);
nRepeats = NaN(nSubj,1);
%shapeStimsRepeated;
k=1; % to iterate thru shapeStimsRepeated

pCorrMandarin = NaN(nSubj,1);
mandarinFam = []; % subjects who speak mandarin
cantoneseFam = []; % subjects who speak cantonese

threshold = NaN(nSubj,1);

pCorrCrm = NaN(nSubj,1);
pCorrMale = NaN(nSubj,1);
pCorrFemale = NaN(nSubj,1);
pCorrCatch = NaN(nSubj,1);

dprime = NaN(nSubj,1);

yrsTrain = NaN(nSubj,1);

for f=1:length(files)
    % load data mat file
    load([dataLoc files{f}]);
    
    % analyze lanna
    lanna = [lanna{1,2},lanna{2,2},lanna{3,2}];
    anchorTrials = lanna(lanna(:,3)==1 | lanna(:,3)==6,1:2);
    testTrials = lanna(lanna(:,3)>1 & lanna(:,3)<6,1:2);
    pCorrAnchor(f) = mean(anchorTrials(:,1)==anchorTrials(:,2));
    pCorrTest(f) = mean(testTrials(:,1)==testTrials(:,2));
    pCorrAll(f) = mean(lanna(:,1)==lanna(:,2));
    
    % analyze Shape
    shapeData = [shape{1,2},shape{2,2},shape{4,2}];
    pCorrShape(f) = mean(shapeData(:,1)==shapeData(:,2));
    nRepeats(f) = sum(shapeData(:,3));
%     if(nRepeats(f)>0)
%         stims = shape{1,3};
%     end
    
    % analyze Mandarin task
    mandarin = [mandarin{1,2},mandarin{2,2}];
    pCorrMandarin(f) = mean(mandarin(:,1)==mandarin(:,2));
    if f>14 && ~isnan(demo{5,2})% first 14 subj weren't asked about familiarity with chinese lang
        for l=1:length(demo{5,2})
            if strcmp(demo{5,2}{l},'Chinese (Mandarin)')
                mandarinFam = [mandarinFam,f];
            end 
            if strcmp(demo{5,2}{l},'Chinese (Cantonese)')
                cantoneseFam = [cantoneseFam,f];
            end
        end
    end
    
    % analyze pd
    pd = [pd{3,2},pd{4,2},pd{1,2},pd{2,2}];
    pd = pd(31:end,:); % last 70 trials
    threshold(f) = PitchCompareWeibull(pd);
    
    % analyze crm
    crm = [crm{1,2},crm{2,2},cell2mat(crm{5,2})];
    maleVoice = crm(~logical(crm(:,3)),1:2);
    pCorrMale(f) = mean(maleVoice(:,1)==maleVoice(:,2));
    femaleVoice = crm(logical(crm(:,3)),1:2);
    pCorrFemale(f) = mean(femaleVoice(:,1)==femaleVoice(:,2));
    pCorrCrm(f) = mean(crm(:,1)==crm(:,2));
    catchTrials = [maleVoice(ismember(maleVoice(:,1),[5;6;19;26],'rows'),1:2);
        femaleVoice(ismember(femaleVoice(:,1),[0;10;18;28],'rows'),1:2)];
    pCorrCatch(f) = mean(catchTrials(:,1)==catchTrials(:,2));
    
    % analyze ts
    type = ts{1,2};
    resp = ts{2,2};
    dprime(f) = analyzeDprime([type(101:150),resp(101:150)]); % last 50 trials
    
    % analyze Golds MSI
    
    % analyze music yrs
    % if snum<15, yrsTrain is idx 5; else idx 6
    yrsTrain(f) = demo{double(f<15)*-1+6,2}; 
    
end

% pCorr lanna data vs dprime
figure; hold on; grid on; box on;
plot(dprime,logit(pCorrAnchor),'o','linewidth',4);
plot(dprime,logit(pCorrTest),'o','linewidth',4);
plot(dprime,logit(pCorrAll),'o','linewidth',4);
% plot([-.8 5],[.5 .5],'k--','linewidth',2)
legend({'Anchor items','Test items','All items'});
title('Lannamaraine Task');
xlabel('Dprime');
ylabel('Logit(Prop corr)');
xlim([-.8 5])
% ylim([0 1])
xTicks = -.5:.5:5;
set(gca,'xtick',xTicks,'fontsize',16,'linewidth',2)

% pCorr shape data vs dprime
figure; grid on; box on; hold on;
plot(dprime,logit(pCorrShape),'o','linewidth',4);
% plot([-.8 5],[.5 .5],'k--','linewidth',2)
title('Shape Task');
xlabel('Dprime');
ylabel('Logit Prop corr');
xlim([-.8 5])
% ylim([0 1])
xTicks = -.5:.5:5;
set(gca,'xtick',xTicks,'fontsize',16,'linewidth',2)

% % nRepeats each subject needed on the Shape task
% figure; grid on; box on; hold on;
% plot(dprime,nRepeats,'o','linewidth',4);
% title('Shape Task, nRepeats');
% xlabel('Dprime');
% ylabel('nRepeats');
% xlim([-.8 5])
% xTicks = -.5:.5:5;
% set(gca,'xtick',xTicks,'fontsize',16,'linewidth',2)

% which shape stims were repeated


% pCorr Mandarin data vs dprime
familiarity_unknown = logical([ones(1,14),zeros(1,length(pCorrMandarin)-14)]);

figure; grid on; box on; hold on;
plot(dprime(familiarity_unknown),logit(pCorrMandarin(familiarity_unknown)),'or','linewidth',4);
plot(dprime(~familiarity_unknown),logit(pCorrMandarin(~familiarity_unknown)),'ob','linewidth',4);
% plot([-.8 5],[.25 .25],'k--','linewidth',2);
title('Mandarin syllables Task');
xlabel('Dprime');
ylabel('Logit Prop corr');
xlim([-.8 5])
% ylim([0 1])
xTicks = -.5:.5:5;
set(gca,'xtick',xTicks,'fontsize',16,'linewidth',2)
legend({'chinese familiarity unknown','no chinese familiarity'});

% pd threshold vs. dprime
plotdpVsThreshold(dprime,threshold);

% pCorr crm data vs dprime
figure; grid on; box on; hold on;
plot(dprime,logit(pCorrCrm),'o','linewidth',4);
plot(dprime,logit(pCorrMale),'o','linewidth',4);
plot(dprime,logit(pCorrFemale),'o','linewidth',4);
plot(dprime,logit(pCorrCatch),'o','linewidth',4);
% plot([-.8 5],[1/32 1/32],'k--','linewidth',2)
legend({'All','Male speaker','Female speaker','Catch Trials'});
title('Multitalker Speech Task');
xlabel('Dprime');
ylabel('Logit Prop corr');
xlim([-.8 5])
% ylim([0 1])
xTicks = -.5:.5:5;
set(gca,'xtick',xTicks,'fontsize',16,'linewidth',2)

plotpCorrvsThreshold(logit(pCorrCrm),threshold)
title('Logit pCorr multitalker task')

plotpCorrvsThreshold(logit(pCorrAll),threshold)
title('Logit pCorr Lannamaraine')

% dprime histogram
figure;
hist(dprime);
xlabel('dprime');
title(sprintf('n = %d',nSubj));

% music training vs dprime
figure;
scatter(yrsTrain,dprime,100,'k','linewidth',4)
axis on
box on
grid on
xlabel('Years of music training');
ylabel('3-task-d^\prime');
xlim([min(yrsTrain)-1 max(yrsTrain)+1])
ylim([-1 5])
% xTicks = [0,5,10,14];
% set(gca,'xtick',xTicks,'fontsize',16,'linewidth',2)
set(gca,'fontsize',16,'linewidth',2)
[r1,p1] = corrcoef(yrsTrain,dprime);
title(sprintf('r = %.2f, p = %.2f',r1(2),p1(2)));
plotRegLine(yrsTrain,dprime,[min(yrsTrain)-1 max(yrsTrain)+1]);

% music training vs threshold
plotMusicVsThreshold(yrsTrain,threshold)

end

function plotdpVsThreshold(TonalityD,thresholdsPC)

thresholdsPC(thresholdsPC<6.25) = 6.25;

figure

log2Thresholds = log2(thresholdsPC);
hold on
plot([-.8 5],-log2(50)*[1 1],'k--','linewidth',2)
scatter(TonalityD,-log2Thresholds,100,'k','linewidth',4)
whichOnes = thresholdsPC > 100 & TonalityD > 1;
plot(TonalityD(whichOnes),-log2Thresholds(whichOnes),'ko','linewidth',4,'markerfacecolor','k')
% yVals = round(2.^(1:.25:3.5));
% yTicks = -log2([1600 800 400 200 100 50 25 12.5 6.25 3.125 1.5625])
yTicks = -log2([1600 800 400 200 100 50 25 12.5 6.25]);
yTickVals = 2.^(-yTicks);
xTicks = -.5:.5:5;
set(gca,'xtick',xTicks,'ytick',yTicks,'yticklabel',yTickVals,'fontsize',16,'linewidth',2)
ylim([-log2(1600) -log2(12.5/2.5)])
xlim([-.8 5])
axis on
box on
grid on
xlabel('3-task-d^\prime')
ylabel('Pitch-difference threshold (cents)')

% edit y-axis label to show <=6.25
labels = strsplit(num2str(yTickVals));
labels{9}='\leq6.25';
yticklabels(labels);

end

function plotMusicVsThreshold(yrsTrain,thresholdsPC)

thresholdsPC(thresholdsPC<6.25) = 6.25;

figure

log2Thresholds = log2(thresholdsPC);
hold on
% plot([-.8 5],-log2(50)*[1 1],'k--','linewidth',2)

scatter(yrsTrain,-log2Thresholds,100,'k','linewidth',4)
% yVals = round(2.^(1:.25:3.5));
% yTicks = -log2([1600 800 400 200 100 50 25 12.5 6.25 3.125 1.5625])
yTicks = -log2([1600 800 400 200 100 50 25 12.5 6.25]);
yTickVals = 2.^(-yTicks);
% xTicks = [0,5,10,14];
% set(gca,'xtick',xTicks,'ytick',yTicks,'yticklabel',yTickVals,'fontsize',16,'linewidth',2)
set(gca,'ytick',yTicks,'yticklabel',yTickVals,'fontsize',16,'linewidth',2)
ylim([-log2(1600) -log2(12.5/2.5)])
axis on
box on
grid on
xlabel('Years of music training');
ylabel('Pitch-difference threshold (cents)')
[r1,p1] = corrcoef(yrsTrain,-log2Thresholds);
plotRegLine(yrsTrain,-log2Thresholds,[min(yrsTrain)-1 max(yrsTrain)+1]);
title(sprintf('r = %.2f, p = %.2f',r1(2),p1(2)));

% edit y-axis label to show <=6.25
labels = strsplit(num2str(yTickVals));
labels{9}='\leq6.25';
yticklabels(labels);

end

function plotpCorrvsThreshold(propCorr,thresholdsPC)

thresholdsPC(thresholdsPC<6.25) = 6.25;

figure

log2Thresholds = log2(thresholdsPC);
hold on
% plot([-.8 5],-log2(50)*[1 1],'k--','linewidth',2)

scatter(propCorr,-log2Thresholds,100,'k','linewidth',4)
% yVals = round(2.^(1:.25:3.5));
% yTicks = -log2([1600 800 400 200 100 50 25 12.5 6.25 3.125 1.5625])
yTicks = -log2([1600 800 400 200 100 50 25 12.5 6.25]);
yTickVals = 2.^(-yTicks);
% xTicks = [0,5,10,14];
% set(gca,'xtick',xTicks,'ytick',yTicks,'yticklabel',yTickVals,'fontsize',16,'linewidth',2)
set(gca,'ytick',yTicks,'yticklabel',yTickVals,'fontsize',16,'linewidth',2)
ylim([-log2(1600) -log2(12.5/2.5)])
axis on
box on
grid on
xlabel('prop correct');
ylabel('Pitch-difference threshold (cents)')
% [r1,p1] = corrcoef(yrsTrain,-log2Thresholds);
% plotRegLine(yrsTrain,-log2Thresholds,[min(yrsTrain)-1 max(yrsTrain)+1]);
% title(sprintf('r = %.2f, p = %.2f',r1(2),p1(2)));

% edit y-axis label to show <=6.25
labels = strsplit(num2str(yTickVals));
labels{9}='\leq6.25';
yticklabels(labels);

end

function plotRegLine(x,y,xlims,color)
if nargin<4
    color = 'k';
end

M = [ones(length(y),1) x(:)];
weights = pinv(M)*y(:);
hold on
plot(xlims,weights(1)+weights(2)*xlims,color,'linewidth',2)
xlim(xlims);
end 