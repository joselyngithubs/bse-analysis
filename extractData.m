function extractData

sourceLoc = 'C:\Users\Joselyn\Documents\GitHub\BSE-data\csv\';
destLoc = 'C:\Users\Joselyn\Documents\GitHub\BSE-data\mat\';

load inits inits
subjNum = length(inits)+1;

files = dir([sourceLoc,'*.csv']);
files = {files.name};

for f=1:length(files)
    
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
    
    % extract Shape data to variable shape
    shapeData = num(num(:,2)==2,4:5);
    shapeStims = txt(num(:,2)==2,4:7);
    shape = {
        'type',shapeData(:,1);
        'resp',shapeData(:,2);
        'stims',shapeStims;
        };
    
    % save initials to inits file
    inits{subjNum} = demo{1,2};
    save inits inits;
    
    % save to mat file
    fname = sprintf('bse%03d',subjNum);
    cd(destLoc);
    if exist([destLoc fname '.mat'])==2
        error('error -- mat file already exists with this subj num');
    else
        save(sprintf('%s',fname),'demo','lanna','shape');
    end
    cd ..
    
    subjNum = subjNum+1;
end

end