function output =  baselineSubtract(input)
% function to baseline subract the first image in a sequence of image
% frames

baselineFrame = repmat(input(:,:,1),[1 1 size(input,3)]);

output = input - baselineFrame;


end