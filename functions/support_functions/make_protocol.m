function protocol = make_protocol(B, Delta, delta)

% Function to create the Camino-like protocol struct for a
% direction-averaged (or spherical mean) signal using a Linear Tensor
% Encoding acquisition scheme (e.g. PGSE)
%
% INPUT
% - B: vector of b values in um^2/ms
% - Delta: gradient pulse separation in ms
% - delta: gradient pulse duration in ms
%
% OUTPUT
% - protocol: the Camino-like protocol structure containing the acquisition
% parameters to be used to generate signal from compartmental models
%
% Author:
% Dr. Marco Palombo
% Cardiff University Brain Research Imaging Centre (CUBRIC)
% Cardiff University, UK
% 8th December 2021
% Email: palombom@cardiff.ac.uk

bx = 1./sqrt(3);
by = 1./sqrt(3);
bz = 1./sqrt(3);
bdir = cell(length(B),1);

for i=1:length(B)
    
    bdir{i} =unique([bx by bz],'rows');
    
end

nbval = length(B);

GAMMA = 2.675987E8;

protocol = struct;

protocol.pulseseq = 'PGSE';
protocol.testrategy = 'fixed';
protocol.roots_sphere = BesselJ_RootsSphere(100);

protocol.delta = Delta.*1e-3.*ones(size(B));
protocol.smalldel = delta.*1e-3.*ones(size(B));

G = sqrt(B.*1e9./(protocol.delta-protocol.smalldel/3))./(GAMMA.*protocol.smalldel);

k=1;

ndir_tot = 0;

for i=1:nbval
    
    ndir = size(bdir{i},1);
    
    protocol.ndir(i) = ndir;
    
    ndir_tot = ndir_tot + ndir;
    
    for j=1:ndir
        protocol.G(k) = G(i);
        protocol.grad_dirs(k,:) = bdir{i}(j,:);
        k=k+1;
    end
end
protocol.ndir_tot = ndir_tot;
protocol.nbval = nbval;

end