function [] = make_direction_average(outputfolder, data_filename, bval_filename, delta, smalldel)

% Calculate the direction-average of the data 
%
% Author:
% Dr. Marco Palombo
% Cardiff University Brain Research Imaging Centre (CUBRIC)
% Cardiff University, UK
% 8th December 2021
% Email: palombom@cardiff.ac.uk

disp('Computing the Direction Average')

tic

mkdir(outputfolder)

% Load the imaging data
tmp = load_untouch_nii(data_filename);
I = tmp.img;
tmp.img = [];

[sx, sy, sz, ~] = size(I);

% Load bvals and bvecs
bvals = importdata(bval_filename);
bvals = round(bvals/100).*100;

bunique = unique(bvals);

Save = zeros(sx, sy, sz, numel(bunique));
S0mean = nanmean(double(I(:,:,:,bvals<=100)),4);

% Identify b-shells and direction-average per shell
for i=1:numel(bunique)
    if i==1
        
        Save(:,:,:,i) = S0mean./S0mean;
        
    else
        
        Save(:,:,:,i) = nanmean(double(I(:,:,:,bvals==bunique(i))),4)./S0mean;
        
    end
end

% Save the direction-averaged data in NIFTI
tmp.img = Save;
tmp.hdr.dime.dim(5) = size(tmp.img,4);
tmp.hdr.dime.dim(1) = 4;
tmp.hdr.dime.datatype = 16;
tmp.hdr.dime.bitpix = 32;

save_untouch_nii(tmp, [outputfolder '/diravg_signal.nii.gz'])

% Make and save the schemefile associated with the direction-averaged data

protocol = make_protocol(bunique./1E3, delta, smalldel);
ProtocolToScheme(protocol, [outputfolder '/diravg.scheme'])
tt = toc;
disp(['DONE - Direction average computed in ' num2str(round(tt)) ' sec.'])

end