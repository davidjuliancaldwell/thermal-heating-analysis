function [] = convert_infratech_csv()
% DJC 7-13-2017 - Conversion file to process Infratech videos
% 
%parpool 

if (exist('myGetenv', 'file'))
    start = myGetenv('subject_dir');    
    if (isempty(start))
        start = pwd;
    end
else
    start = pwd;
end

[filename,rawpath] = uigetfile(start, 'select a beginning file to convert');
cd(rawpath)
D = dir([rawpath, '\*.csv']);
Num = length(D(not([D.isdir])));

files = dir([rawpath,'\*.csv']);

ind = 1;
%stackedData = zeros(480,640,Num); % this is if using entire screen oof caemra 
% whole image - 480 x 640 pixels in each frame
wholeImage = 1;

for file = files'
    
    data_single = importfile_infra(file.name,wholeImage);
    
    if ind == 1
       stackedData = zeros(size(data_single,1),size(data_single,2),Num); 
    end
    
    stackedData(:,:,ind) = data_single;
    ind = ind+1;
    fprintf(['file ' num2str(ind) '\n'])

end

split_name = strsplit(filename,'.csv');

save(fullfile(rawpath, [split_name{1} '.mat']), '-v7.3', 'stackedData');

end

