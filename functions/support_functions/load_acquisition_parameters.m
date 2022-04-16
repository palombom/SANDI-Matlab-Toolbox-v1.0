function [] = load_acquisition_parameters(bvalues_filename, delta, smalldelta, outputfolder)
% Function to load ad save the acquisition parameters in camino-like
% scheme file
%
% Author:
% Dr. Marco Palombo
% Cardiff University Brain Research Imaging Centre (CUBRIC)
% Cardiff University, UK
% 8th December 2021
% Email: palombom@cardiff.ac.uk

disp('Loading the acquisition parameters')
bvals = importdata(bvalues_filename);
bvals = round(bvals/100).*100; % to round the b value in case it changes slightly across directions
bunique = unique(bvals);
protocol = make_protocol(bunique./1E3, delta, smalldelta);
schemefile = fullfile(outputfolder, 'diravg.scheme');
ProtocolToScheme(protocol, schemefile)

end