function [ annotation_vector ] = get_annotation_vector( max_length, annSamples, annComments )
    annotation_vector = zeros(max_length, 1);

    for i=1:(length(annComments))    
        % Start of the range
        startSample = annSamples(i);

        % End of the range
        if i~=length(annComments)
            endSample = annSamples(i+1)-1;    
        else
            endSample = max_length;
        end

        % Write on annotation_vector
        if strcmp(annComments(i), '(AFIB')
            annotation_vector(startSample:endSample, 1) = 1;
        elseif ~strcmp(annComments(i), '(N')
            annotation_vector(startSample:endSample, 1) = 2;
        end
    end
end

