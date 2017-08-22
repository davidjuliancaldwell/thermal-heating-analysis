function dataOut = importfile_infra(filename)

% e.g. Sequencefile01 = convert_demo('Sequencefile_01.csv');

%% Open the text file.
fileID = fopen(filename,'r');

delimiter = ';';
wholeImage = 1;
if wholeImage
    startRow = 9;
    
else
    startRow = 11;
end
textscan(fileID, '%[^\n\r]', startRow(1)-1, 'WhiteSpace', '', 'ReturnOnError', false);
dataArray = textscan(fileID, '%[^\n\r]','Delimiter', delimiter, 'ReturnOnError', false);

fclose(fileID);

dataArray = cell2mat(dataArray{1});
%dataOut = zeros(size(dataArray,1),640);
%dataOut = zeros(size(dataArray,1),size(dataArray,2)/4);
    i = 1;
    temp = dataArray(i,:);
    temp_split = strsplit(temp,delimiter);
    temp_num_cell = strrep(temp_split,',','.');
    dataOut = zeros(size(dataArray,1),size(temp_num_cell,2)-1); 
    dataOut(i,:) = str2double(temp_num_cell(:,1:end-1));

parfor i = 2:size(dataArray,1)
    
    temp = dataArray(i,:);
    temp_split = strsplit(temp,delimiter);
    temp_num_cell = strrep(temp_split,',','.');

    dataOut(i,:) = str2double(temp_num_cell(:,1:end-1));
end
end