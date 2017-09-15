function output =  baselineSubtract(input)
% DJC BLBT Summer 2017. This is a function to subtract the first image in a
% whole sequence from all sequential frames for further analysis in terms
% of the heating videos

baselineFrame = repmat(input(:,:,1),[1 1 size(input,3)]);

output = input - baselineFrame;


end