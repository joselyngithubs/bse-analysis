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
shape_stims_repeated = cell(1,1);
k=1; % to iterate thru shape_stims_repeated
idx_set = {'A','B','C'};

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

gmsiScores = NaN(nSubj,3);

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
    % collect which stimuli were repeated
    if(nRepeats(f)>0)
        stims_tmp = shape{3,2};
        n_repeated = shapeData(:,3);
        
        % filter just the stims that were repeated
        stims_tmp = stims_tmp(n_repeated>0,:);
        n_repeated = n_repeated(n_repeated>0);
        
        for r=1:size(stims_tmp,1)
            % the last column is A,B,C so convert those to numbers first
            idx_tmp = find(strcmp(stims_tmp(r,4),idx_set));
            % grab the relevant word based on that index; store in
            % shape_stims_repeated as many times as that word was repeated
            for m=1:n_repeated(r)
                shape_stims_repeated{k} = stims_tmp{idx_tmp};
                k = k+1;
            end
            
        end
    end
    
    % analyze Mandarin task
    mandarin = [mandarin{1,2},mandarin{2,2}];
    pCorrMandarin(f) = mean(mandarin(:,1)==mandarin(:,2));
    if f>14 % first 14 subj weren't asked about familiarity with chinese lang
        if length(demo{5,2})>1 % for subj 15+, check if demo{5,2} contains something
        % if demo{5,2} contains something, it will either be "Chinese
        % (Cantonese)", "Chinese (Mandarin)", or BOTH of those strings as
        % one combined string. So, check the contents based on length of
        % the full string.
            if length(demo{5,2})>19 % if longer than "Chinese (Cantonese)", then both languages must be listed
                cantoneseFam = [cantoneseFam,f];           
                mandarinFam = [mandarinFam,f];
            elseif length(demo{5,2})<19 % only contains "Chinese (Mandarin)"
                mandarinFam = [mandarinFam,f];
            else % only contains "Chinese (Cantonese)"
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
    gmsiScores(f,:) = [gmsi{4,2}(10),gmsi{4,2}(18),gmsi{4,2}(37)];
    
    % analyze music yrs
    % if snum<15, yrsTrain is idx 5; else idx 6
    yrsTrain(f) = demo{double(f<15)*-1+6,2}; 
    
end

%% for the dissertation
% makeDissertationPlot(dprime,threshold,pCorrShape);

%% plot Lanna stat
lannaStat = getLannaStat;

thresholdsPC = threshold;
thresholdsPC(thresholdsPC<6.25) = 6.25;
log2Thresholds = log2(thresholdsPC);
figure;
hold on
scatter(lannaStat,-log2Thresholds,100,'k','linewidth',4)
yTicks = -log2([1600 800 400 200 100 50 25 12.5 6.25]);
yTickVals = 2.^(-yTicks);
set(gca,'ytick',yTicks,'yticklabel',yTickVals,'fontsize',16,'linewidth',2)
ylim([-log2(1600) -log2(12.5/2.5)])
axis on
box on
grid on
xlabel('lanna stat');
ylabel('Pitch-difference threshold (cents)')
% edit y-axis label to show <=6.25
labels = strsplit(num2str(yTickVals));
labels{9}='\leq6.25';
yticklabels(labels);
[r1,p1] = corrcoef(lannaStat,-log2Thresholds);
title(sprintf('r = %.2f, p = %.2f',r1(2),p1(2)));
plot([0 0],[-log2(1600) -log2(12.5/2.5)],'k--','linewidth',2)

% figure;
% histogram(logit(pCorrAll),7);
% title('lanna logit pcorr')
% 
% figure;
% histogram(lannaStat,7);
% title('lanna stat');


%%
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

% % pCorr shape data vs dprime
% figure; grid on; box on; hold on;
% plot(dprime,logit(pCorrShape),'o','linewidth',4);
% % plot([-.8 5],[.5 .5],'k--','linewidth',2)
% title('Shape Task');
% xlabel('Dprime');
% ylabel('Logit Prop corr');
% xlim([-.8 5])
% % ylim([0 1])
% xTicks = -.5:.5:5;
% set(gca,'xtick',xTicks,'fontsize',16,'linewidth',2)

% pCorr shape data vs dprime (legend to show if threshold below or above 50
% cents)
figure; grid on; box on; hold on;
whichOnes = threshold < 50;
plot(dprime(whichOnes),logit(pCorrShape(whichOnes)),'ko','linewidth',4);
plot(dprime(~whichOnes),logit(pCorrShape(~whichOnes)),'ro','linewidth',4);
% plot([-.8 5],[.5 .5],'k--','linewidth',2)
xlabel('3-task-d^\prime');
ylabel('Logit Prop corr');
xlim([-.8 5])
% ylim([0 1])
xTicks = -.5:.5:5;
set(gca,'xtick',xTicks,'fontsize',16,'linewidth',2)
legend({'PDT < 50','PDT > 50'});
[r1,p1] = corrcoef(dprime(whichOnes),logit(pCorrShape(whichOnes)));
[r2,p2] = corrcoef(dprime(~whichOnes),logit(pCorrShape(~whichOnes)));
title('Shape Task')
text(2,-1.1,sprintf('below50: r = %.2f, p = %.2f\nabove50: r = %.2f, p = %.2f',r1(2),p1(2),r2(2),p2(2)));

% % for the defense slides
% % pCorr shape data vs dprime (legend to show if threshold below or above 50
% % cents)
% gold_color = [218,165,32];
% figure; grid on; box on; hold on;
% whichOnes = threshold < 50;
% plot(dprime(whichOnes),logit(pCorrShape(whichOnes)),'ko','linewidth',4);
% plot(dprime(~whichOnes),logit(pCorrShape(~whichOnes)),'o','linewidth',4,'markeredgecolor',gold_color./218);
% % plot([-.8 5],[.5 .5],'k--','linewidth',2)
% xlabel('3-task-d^\prime');
% ylabel('Prosody score');
% xlim([-.8 5])
% % ylim([0 1])
% xTicks = -.5:.5:5;
% set(gca,'xtick',xTicks,'fontsize',16,'linewidth',2)
% legend({'PDT < 50','PDT > 50'});
% [r1,p1] = corrcoef(dprime(whichOnes),logit(pCorrShape(whichOnes)));
% [r2,p2] = corrcoef(dprime(~whichOnes),logit(pCorrShape(~whichOnes)));

% % nRepeats each subject needed on the Shape task
figure; grid on; box on; hold on;
plot(dprime,nRepeats,'o','linewidth',4);
title('Shape Task, nRepeats');
xlabel('Dprime');
ylabel('nRepeats');
xlim([-.8 5])
xTicks = -.5:.5:5;
set(gca,'xtick',xTicks,'fontsize',16,'linewidth',2)

% which shape stims were repeated
figure
histogram(categorical(shape_stims_repeated));
title('number of times each Shape stim was repeated');
xtickangle(45);

%% Mandarin
% pCorr Mandarin data vs dprime. Note that there may be data plotted twice
% if a subject is familiar with BOTH cantonese and mandarin.
familiarity_unknown = logical([ones(1,14),zeros(1,length(pCorrMandarin)-14)]);
no_chinese = ~familiarity_unknown;
no_chinese(cantoneseFam) = zeros(size(cantoneseFam));
no_chinese(mandarinFam) = zeros(size(mandarinFam));

figure; grid on; box on; hold on;
plot(dprime(familiarity_unknown),logit(pCorrMandarin(familiarity_unknown)),'ok','linewidth',4);
plot(dprime(no_chinese),logit(pCorrMandarin(no_chinese)),'ob','linewidth',4);
plot(dprime(cantoneseFam),logit(pCorrMandarin(cantoneseFam)),'or','linewidth',4);
plot(dprime(mandarinFam),logit(pCorrMandarin(mandarinFam)),'og','linewidth',4);
% plot([-.8 5],[.25 .25],'k--','linewidth',2);
title('Mandarin syllables Task');
xlabel('Dprime');
ylabel('Logit Prop corr');
xlim([-.8 5])
% ylim([0 1])
xTicks = -.5:.5:5;
set(gca,'xtick',xTicks,'fontsize',16,'linewidth',2)
legend({'chinese familiarity unknown','no chinese familiarity','cantonese familiarity','mandarin familiarity'});

% pCorr Mandarin vs threshold, with color legend. Note that there may be data plotted twice
% if a subject is familiar with BOTH cantonese and mandarin.
thresholdsPC = threshold;
thresholdsPC(thresholdsPC<6.25) = 6.25;
log2Thresholds = log2(thresholdsPC);
figure;
hold on
scatter(logit(pCorrMandarin(familiarity_unknown)),-log2Thresholds(familiarity_unknown),100,'k','linewidth',4)
scatter(logit(pCorrMandarin(no_chinese)),-log2Thresholds(no_chinese),100,'b','linewidth',4)
scatter(logit(pCorrMandarin(cantoneseFam)),-log2Thresholds(cantoneseFam),100,'r','linewidth',4)
scatter(logit(pCorrMandarin(mandarinFam)),-log2Thresholds(mandarinFam),100,'g','linewidth',4)
legend({'chinese familiarity unknown','no chinese familiarity','cantonese familiarity','mandarin familiarity'});
yTicks = -log2([1600 800 400 200 100 50 25 12.5 6.25]);
yTickVals = 2.^(-yTicks);
set(gca,'ytick',yTicks,'yticklabel',yTickVals,'fontsize',16,'linewidth',2)
ylim([-log2(1600) -log2(12.5/2.5)])
axis on
box on
grid on
xlabel('logit prop correct');
ylabel('Pitch-difference threshold (cents)')
% edit y-axis label to show <=6.25
labels = strsplit(num2str(yTickVals));
labels{9}='\leq6.25';
yticklabels(labels);
% line at chance performance 
plot([logit(1/4) logit(1/4)],[-log2(1600) -log2(12.5/2.5)],'k--','linewidth',2)

%%

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

plotpCorrvsThreshold(pCorrCrm,threshold,1/32)
title('Logit pCorr multitalker task')

% plotpCorrvsThreshold(pCorrAll,threshold,1/6)
% title('Logit pCorr Lannamaraine')

plotpCorrvsThreshold(pCorrShape,threshold,1/3)
title('Logit pCorr Shape')

% commenting this out bc we already have the color legend version of this
% to separate people with chinese familiarity
% plotpCorrvsThreshold(pCorrMandarin,threshold,1/4)
% title('Logit pCorr Mandarin')

%% correlation between individual tasks (all logit)

figure;
subplot(1,3,1)
correlateTasks(logit(pCorrShape),logit(pCorrAll));
xlabel('shape');ylabel('lannamaraine')

subplot(1,3,2)
correlateTasks(logit(pCorrMandarin),logit(pCorrAll));
xlabel('mandarin');ylabel('lannamaraine')

subplot(1,3,3)
correlateTasks(logit(pCorrShape),logit(pCorrMandarin));
xlabel('shape');ylabel('mandarin')


%% dprime histogram
% figure;
% hist(dprime);
% xlabel('dprime');
% title(sprintf('n = %d',nSubj));

% plot hist of dp along with what the dp distribution would look like if
% remove thresholds above 50 cents
TonalityD = dprime;
thresholdsPC = threshold;
figure
whichOnes = thresholdsPC < inf;
myhist=histogram(TonalityD(whichOnes),10,'linewidth',2,'facecolor',[.7 .7 .7],'facealpha',1);
get(myhist)
whichOnes = thresholdsPC < 50;
hold on
histogram(TonalityD(whichOnes),'binedges',myhist.BinEdges+.04,'facecolor',[1 1 1],'facealpha',.7,'linewidth',2)
xtick=-1:1:4;
set(gca,'linewidth',2,'fontsize',18,'ytick',0:5:40)
grid on
box on
xlim([-1,4.5])
ylim([0,35])
xticks(xtick)
legend({'All listeners',sprintf('Listeners with\nthreshold <50cents')});
ylabel('Number of listeners')
xlabel('3-task-d^\prime')



%% music training

% music training vs dprime
figure;
subplot(1,2,1)
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
subplot(1,2,2)
plotMusicVsThreshold(yrsTrain,threshold)

%% GMSI
% col 1: perceptual subscale score
% col 2: training subscale score
% col 3: sophistication subscale score

% gmsi vs threshold
% figure;
% plotGMSIvsThreshold(gmsiScores(:,1),threshold);
% title('GMSI perceptual score')
% figure;
% plotGMSIvsThreshold(gmsiScores(:,2),threshold);
% title('GMSI training score')
figure;
subplot(1,3,1)
plotGMSIvsThreshold(gmsiScores(:,3),threshold);
[r1,p1] = corrcoef(gmsiScores(:,3),threshold);
title(sprintf('r = %.2f, p = %.2f',r1(2),p1(2)));

% GMSI sophistication vs dprime
subplot(1,3,2)
scatter(gmsiScores(:,3),dprime,100,'k','linewidth',4)
axis on
box on
grid on
xlabel('GMSI sophistication score');
ylabel('3-task-d^\prime');
xlim([min(gmsiScores(:,3))-1 max(gmsiScores(:,3))+1])
ylim([-1 5])
% xTicks = [0,5,10,14];
% set(gca,'xtick',xTicks,'fontsize',16,'linewidth',2)
set(gca,'fontsize',16,'linewidth',2)
[r1,p1] = corrcoef(gmsiScores(:,3),dprime);
title(sprintf('r = %.2f, p = %.2f',r1(2),p1(2)));
plotRegLine(gmsiScores(:,3),dprime,[min(gmsiScores(:,3))-1 max(gmsiScores(:,3))+1]);

% % GMSI perceptual vs dprime
% figure;
% scatter(gmsiScores(:,1),dprime,100,'k','linewidth',4)
% axis on
% box on
% grid on
% xlabel('GMSI perceptual score');
% ylabel('3-task-d^\prime');
% xlim([min(gmsiScores(:,1))-1 max(gmsiScores(:,1))+1])
% ylim([-1 5])
% % xTicks = [0,5,10,14];
% % set(gca,'xtick',xTicks,'fontsize',16,'linewidth',2)
% set(gca,'fontsize',16,'linewidth',2)
% [r1,p1] = corrcoef(gmsiScores(:,1),dprime);
% title(sprintf('r = %.2f, p = %.2f',r1(2),p1(2)));
% plotRegLine(gmsiScores(:,1),dprime,[min(gmsiScores(:,1))-1 max(gmsiScores(:,1))+1]);

% % GMSI training vs dprime
% figure;
% scatter(gmsiScores(:,2),dprime,100,'k','linewidth',4)
% axis on
% box on
% grid on
% xlabel('GMSI training score');
% ylabel('3-task-d^\prime');
% xlim([min(gmsiScores(:,2))-1 max(gmsiScores(:,2))+1])
% ylim([-1 5])
% % xTicks = [0,5,10,14];
% % set(gca,'xtick',xTicks,'fontsize',16,'linewidth',2)
% set(gca,'fontsize',16,'linewidth',2)
% [r1,p1] = corrcoef(gmsiScores(:,2),dprime);
% title(sprintf('r = %.2f, p = %.2f',r1(2),p1(2)));
% plotRegLine(gmsiScores(:,2),dprime,[min(gmsiScores(:,2))-1 max(gmsiScores(:,2))+1]);

% % gmsi vs gmsi (correlations)
% figure;plot(gmsiScores(:,1),gmsiScores(:,2),'o');
% xlabel('GMSI perceptual');ylabel('GMSI training');
% [r1,p1] = corrcoef(gmsiScores(:,1),gmsiScores(:,2));
% title(sprintf('r = %.2f, p = %.2f',r1(2),p1(2)));
% plotRegLine(gmsiScores(:,1),gmsiScores(:,2),[min(gmsiScores(:,1))-1 max(gmsiScores(:,1))+1]);
% figure;plot(gmsiScores(:,1),gmsiScores(:,3),'o');
% xlabel('GMSI perceptual');ylabel('GMSI sophistication');
% [r1,p1] = corrcoef(gmsiScores(:,1),gmsiScores(:,3));
% title(sprintf('r = %.2f, p = %.2f',r1(2),p1(2)));
% plotRegLine(gmsiScores(:,1),gmsiScores(:,3),[min(gmsiScores(:,1))-1 max(gmsiScores(:,1))+1]);
% figure;plot(gmsiScores(:,3),gmsiScores(:,2),'o');
% xlabel('GMSI sophistication');ylabel('GMSI training');
% [r1,p1] = corrcoef(gmsiScores(:,3),gmsiScores(:,2));
% title(sprintf('r = %.2f, p = %.2f',r1(2),p1(2)));
% plotRegLine(gmsiScores(:,3),gmsiScores(:,2),[min(gmsiScores(:,3))-1 max(gmsiScores(:,3))+1]);

subplot(1,3,3)
plot(gmsiScores(:,3),yrsTrain,'o');
ylabel('Yrs music training')
xlabel('GMSI sophistication score')
[r1,p1] = corrcoef(gmsiScores(:,3),yrsTrain);
title(sprintf('r = %.2f, p = %.2f',r1(2),p1(2)));


end

function plotGMSIvsThreshold(gmsi,thresholdsPC)

thresholdsPC(thresholdsPC<6.25) = 6.25;

log2Thresholds = log2(thresholdsPC);
hold on
% plot([-.8 5],-log2(50)*[1 1],'k--','linewidth',2)
scatter(gmsi,-log2Thresholds,100,'k','linewidth',4)

yTicks = -log2([1600 800 400 200 100 50 25 12.5 6.25]);
yTickVals = 2.^(-yTicks);
% xTicks = -.5:.5:5;
set(gca,'ytick',yTicks,'yticklabel',yTickVals,'fontsize',16,'linewidth',2)
ylim([-log2(1600) -log2(12.5/2.5)])
% xlim([-.8 5])
axis on
box on
grid on
xlabel('gmsi score')
ylabel('Pitch-difference threshold (cents)')

% edit y-axis label to show <=6.25
labels = strsplit(num2str(yTickVals));
labels{9}='\leq6.25';
yticklabels(labels);

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

function plotpCorrvsThreshold(propCorr,thresholdsPC,chancePerformance)

propCorr = logit(propCorr);

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
xlabel('logit prop correct');
ylabel('Pitch-difference threshold (cents)')
% [r1,p1] = corrcoef(yrsTrain,-log2Thresholds);
% plotRegLine(yrsTrain,-log2Thresholds,[min(yrsTrain)-1 max(yrsTrain)+1]);
% title(sprintf('r = %.2f, p = %.2f',r1(2),p1(2)));

% line at chance performance depending on the task
plot([logit(chancePerformance) logit(chancePerformance)],[-log2(1600) -log2(12.5/2.5)],'k--','linewidth',2)

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

function makeDissertationPlot(TonalityD,thresholdsPC,propCorr)

propCorr = logit(propCorr);
thresholdsPC(thresholdsPC<6.25) = 6.25;

figure('Renderer', 'painters', 'units','normalized','Position', [0 .3 .75 .6])
y = .2;
wd = .4;
ht = .75;
subplot('Position',[.12 y wd ht]);

log2Thresholds = log2(thresholdsPC);
hold on
plot([-.8 5],-log2(50)*[1 1],'k--','linewidth',2)
scatter(TonalityD,-log2Thresholds,100,'k','linewidth',4)
whichOnes = thresholdsPC <0;
% whichOnes = thresholdsPC > 100 & TonalityD > 1;
plot(TonalityD(whichOnes),-log2Thresholds(whichOnes),'ko','linewidth',4,'markerfacecolor','k')
% yVals = round(2.^(1:.25:3.5));
% yTicks = -log2([1600 800 400 200 100 50 25 12.5 6.25 3.125 1.5625])
yTicks = -log2([1600 800 400 200 100 50 25 12.5 6.25]);
yTickVals = 2.^(-yTicks);
xTicks = -.5:.5:5;
xTickVals = {'','0','','1','','2','','3','','4','','5'};
set(gca,'xtick',xTicks,'ytick',yTicks,'yticklabel',yTickVals,'xticklabel',xTickVals,'fontsize',16,'linewidth',2)
ylim([-log2(1600) -log2(12.5/2.5)])
xlim([-.8 5])
axis on
box on
grid on
xlabel('3-task-d^\prime')
ylabel('PDT (cents)')

% edit y-axis label to show <=6.25
labels = strsplit(num2str(yTickVals));
labels{9}='\leq6.25';
yticklabels(labels);

subplot('Position',[.14+wd y wd ht]);

scatter(propCorr,-log2Thresholds,100,'k','linewidth',4)
hold on
yTicks = -log2([1600 800 400 200 100 50 25 12.5 6.25]);
yTickVals = 2.^(-yTicks);
xTicks = -2:1:4;
set(gca,'xtick',xTicks,'ytick',yTicks,'yticklabel',{},'fontsize',16,'linewidth',2)
ylim([-log2(1600) -log2(12.5/2.5)])
axis on
box on
grid on
xlabel('Shape Task Score');
plot([-2 4],-log2(50)*[1 1],'k--','linewidth',2)

end

function correlateTasks(x,y)
% note that this function is tailored for comparison of Lanna, Mandarin,
% and Shape tasks only.

plot(x,y,'o');
ylim([-3,5])
xlim([-3,5])

% plot reg line
[r1,p1] = corrcoef(x,y);
plotRegLine(x,y,[-3,5]);
title(sprintf('r = %.2f, p = %.2f',r1(2),p1(2)));

end