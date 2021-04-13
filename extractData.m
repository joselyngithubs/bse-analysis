function extractData

sourceLoc = 'C:\Users\Joselyn\Documents\GitHub\BSE-data\csv\';
destLoc = 'C:\Users\Joselyn\Documents\GitHub\BSE-data\mat\';

load inits inits
subjNum = length(inits)+1;

files = dir([sourceLoc,'*.csv']);
files = {files.name};

for f=1:length(files)
    fprintf('\nSubject %d of %d',f,length(files));
    
    % extract demographic data to variable demo
    [~,~,demo] = xlsread([sourceLoc,files{f}],'A2:H2');
    [taskOrder,~,~] = xlsread([sourceLoc,files{f}],'J2:J544');
    taskOrder = unique(taskOrder,'stable');
    demo = {
        'snum',subjNum;
        'init',demo{1};
        'lang',demo{3};
        'langOther',demo{4};
        'chinese',demo{5};
        'yrsTrain',demo{6};
        'fs',demo{7};
        'taskOrder',taskOrder
        };
    
    % extract all task data
    [num,txt,raw] = xlsread([sourceLoc,files{f}],'I2:R544');
    
    % extract Lanna data (taskID 1) to variable lanna
    lannaData = num(num(:,2)==1,[4,5,6]);
    lanna = {
        'type',lannaData(:,1);
        'resp',lannaData(:,2);
        'stim',lannaData(:,3);
        };
    
    % extract Shape data (taskID 2) to variable shape
    shapeData = num(num(:,2)==2,[4:5,10]);
    shapeStims = txt(num(:,2)==2,4:7);
    shape = {
        'type',shapeData(:,1);
        'resp',shapeData(:,2);
        'stims',shapeStims;
        'nRepeat',shapeData(:,3);
        };
    
    % extract Mandarin data (taskID 3) to variable mandarin
    mandarinData = num(num(:,2)==3,4:5);
    mandarinStims = txt(num(:,2)==3,4);
    mandarin = {
        'type',mandarinData(:,1);
        'resp',mandarinData(:,2);
        'stims',mandarinStims;
        };
    
    % extract pitch difference data (taskID 4) to variable pd
    pdData = num(num(:,2)==4,4:7);
    pd = {
        'type',pdData(:,1);
        'resp',pdData(:,2);
        'freq1',pdData(:,3);
        'freq2',pdData(:,4);
        };
    
    % extract CRM data (taskID 5) to variable crm
    crmData = num(num(:,2)==5,4:5);
    crmStims = raw(num(:,2)==5,6:8);
    crm = {
        'type',crmData(:,1);
        'resp',crmData(:,2);
        'stims',crmStims(:,1);
        'stimIndex',crmStims(:,2);
        'speaker',crmStims(:,3);
    };

    % extract tone-scramble data (taskID 6) to variable ts
    tsData = num(num(:,2)==6,4:5);
    tsStims = raw(num(:,2)==6,6);
    ts = {
        'type',tsData(:,1);
        'resp',tsData(:,2);
        'seq',tsStims(:,1);
        };
    
    % extract Golds MSI data (taskID 7) to variable gmsi
    gmsiData = num(num(:,2)==7,[4:6,8]);
    gmsiStims = txt(num(:,2)==7,5);
    gmsi = {
        'type',gmsiData(:,1);
        'resp',gmsiData(:,2);
        'index',gmsiData(:,3);
        'score',gmsiData(:,4);
        'subscale',gmsiStims;
    };
    
    % save initials to inits file
    inits{subjNum} = demo{2,2};
    save inits inits;
    
    % save to mat file
    fname = sprintf('bse%03d',subjNum);
    cd(destLoc);
    if exist([destLoc fname '.mat'])==2
        error('error -- mat file already exists with this subj num');
    else
        save(sprintf('%s',fname),'demo','lanna','shape','mandarin','pd','crm','ts','gmsi');
    end
    cd ..
    
    subjNum = subjNum+1;
end

end