function mpgMean = apply_GRNN_matlab(signal, Mdl)

% Apply the Generalized Regression Neural Network
%
% Author:
% Dr. Marco Palombo
% Cardiff University Brain Research Imaging Centre (CUBRIC)
% Cardiff University, UK
% 8th December 2021
% Email: palombom@cardiff.ac.uk

tic

input_size = size(signal,1);
batch_size = 7000;

if input_size>batch_size
    
    tmp = 1:batch_size:input_size;
    n_batches = numel(tmp);
    
    mpgMean = zeros(5, input_size);
    
    for i = 1:n_batches-1
        
        mpgMean(:, tmp(i):tmp(i)+batch_size-1) = Mdl(signal(tmp(i):tmp(i)+batch_size-1, :)');
        
    end
    
    mpgMean(:, tmp(end):end) = Mdl(signal(tmp(end):end, :)');
    
    
else
    
    mpgMean = Mdl(signal');
    
end

tt = toc;

mpgMean = mpgMean';

fprintf('DONE - RF fitted in %3.0f sec.\n', tt)

end