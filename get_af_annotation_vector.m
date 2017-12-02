function [ annotation_vector ] = get_af_annotation_vector( max_length, annSamples, annComments )
%GET_ANNOTATION_VECTOR Returns a single Nx1 logical vector containing 1
%for samples marked as being part of an atrial fibrillation sequence

annotation_vector = false(max_length, 1);

% Handle all comments except the last
for i=1:(length(annComments)-1)    
    % We only care about (AFIB comments
    if strcmp(annComments(i), '(AFIB')
        startSample = annSamples(i);
        endSample = annSamples(i+1)-1;

        annotation_vector(startSample:endSample, 1) = 1;
    end
end

% Handle situation where the very last comment is (AFIB
if(strcmp(annComments(end), '(AFIB')) 
    startSample = annSamples(end);
    endSample = max_length;
    annotation_vector(startSample:endSample, 1) = 1;
end

end

