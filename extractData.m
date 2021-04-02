sourceLoc = 'C:\Users\Joselyn\Documents\GitHub\BSE-data\csv\';
destLoc = 'C:\Users\Joselyn\Documents\GitHub\BSE-data\mat';

cd(destLoc)
subjects = dir('*.mat'); % count # of subjects to assign subj number
subjNum = length(subjects)+1;
cd ..

files = dir([sourceLoc,'*.csv']);
files = {files.name};

demo = cell(1,7);
lanna = NaN(96,10);

for f=1:length(files)
    fname = sprintf('bse%03d',subjNum);
    subjNum = subjNum+1;
    
    % extract demographic data to variable demo
    [~,~,demo] = xlsread([sourceLoc,files{f}],'A2:G2');
    [taskOrder,~,~] = xlsread([sourceLoc,files{f}],'I2:I544');
    taskOrder = unique(taskOrder,'stable');
    demo = {
        'init',demo{1};
        'lang',demo{3};
        'langOther',demo{4};
        'yrsTrain',demo{5};
        'fs',demo{6};
        'taskOrder',taskOrder
        };
    
    % extract all task data
    [num,txt,raw] = xlsread([sourceLoc,files{f}],'H2:Q544');
    
    % extract Lanna data to variable lanna
    lannaData = num(num(:,2)==1,[4,5,6]);
    lanna = {
        'type',lannaData(:,1);
        'resp',lannaData(:,2);
        'stim',lannaData(:,3);
        };
    
    % save to mat file
    % first check that it doesnt already exist
    save fname demo lanna % this saves it literally as "fname.mat"
end