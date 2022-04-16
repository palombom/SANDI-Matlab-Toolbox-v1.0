function [E,J]=AstroSticks(x, protocol,x_deriv)
%
% camino.m--------------------------------------------------------------
% Diffusion signal simulating the AstroSticks compartment.
% 
% [E,J]=AstroSticks(x, protocol,x_deriv)
% 
% Description: Returns the diffusion measurement and the Jacobian with 
% respect to model parameters for a substrate consisting of isotropically 
% oriented sticks and a diffusion protocol specified in the input
% Substrate: Isotropically distributed sticks
% Diffusion pulse sequence: wrapper for various sequences
%
% Parameters:   
% E - diffusion signal
% J - Jacobian of the diffusion signal with respect to model parameters
% x - size 1 vector of model parameters in SI units for AstroSticks:
%       x(1) - free diffusivity of the material along the sticks     
% protocol - structure which includes all the information related to the 
%        diffusion protocol and required to generate the signal.
%       The basic information for a PGSE sequence is the following:
%       protocol.grad_dirs - is the gradient direction for each measurement.
%           It has size [N 3] where N is the number of measurements.
%       protocol.G - gradient strength, size [1 N]
%       protocol.delta - pulse separation, size [1 N]
%       protocol.smalldel - pulse duration, size [1 N]
%       Other diffusion sequences might have additional fields.
% x_deriv - a vector of 0s and 1s with the same size as x, indicating
%       which parameters are considered in the Jacobian;
%
%------------------------------------------------------------------------
% This file is part of the camino.m toolbox.
% Copyright (c) 2015, UCL Microstructure Imaging Group (MIG), All rights reserved.
% Distributed under the Modified BSD Licence (see: LICENSE.pdf).
%
% Authors:
%   Andrada Ianus (a.ianus.11@ucl.ac.uk)
%   Daniel C. Alexander (d.alexander@ucl.ac.uk)
%      

if strcmp(protocol.pulseseq,'PGSE') || strcmp(protocol.pulseseq,'SWOGSE') || strcmp(protocol.pulseseq,'TWOGSE') || ...
        strcmp(protocol.pulseseq,'DSE') || strcmp(protocol.pulseseq,'STEAM')% sequences with one gradient orientation

    fun = str2func('Stick');
    G = protocol.G';

    protocol1 = protocol; %  make protocol with gradient gradient oriented along z direction
    protocol1.grad_dirs = repmat([0 0 1],length(protocol.smalldel),1);

    Lperp = zeros(size(G)); % no perpendicular component
    Lpar = log(fun([x 0 0],protocol1)); % compute parallel component (cylinder along z, gradient along z)
    
    E_r = ones(size(G));
    E_r(G>0) = sqrt(pi)./(2.*sqrt(Lperp(G>0)-Lpar(G>0))).*exp(Lperp(G>0)).*erf(sqrt(Lperp(G>0)-Lpar(G>0)));    

    E=E_r;

    % Compute the Jacobian matrix; computed numerically
    if(nargout>1)
        dx = 0.00001;
        J = zeros(length(E), 1);
        if nargin < 3          
            for i = 1; % compute the derivatives for all model parameters
                xpert = x;
                xpert(i) = xpert(i)*(1+dx);
                Epert = AstroSticks(xpert, protocol);
                dEdx = (Epert - E)/(xpert(i)*dx);
                J(:,i) = dEdx;
            end
        else % x_deriv is a vector containing 1 and 0  which specifies the parameters that enter the jacobian

            for i = 1:length(x_deriv);
                if x_deriv(i) ~= 0                  
                    xpert = x;
                    xpert(i) = xpert(i)*(1+dx);
                    Epert = AstroSticks(xpert, protocol);
                    dEdx = (Epert - E)/(xpert(i)*dx);
                    J(:,i) = dEdx;
                end
            end
        end   
    end
else % sequences that have gradients with multiple orientations
    % need to loop over different directions
    % get Ndir directions on a spehere
    fun = str2func('Stick');
    Ndir = 100;    
    fibre_dirs = ReadCaminoElecPS(sprintf('PointSets/Elec%03i.txt',Ndir))';  
    [az incl] = cart2sph(fibre_dirs(1,:),fibre_dirs(2,:),fibre_dirs(3,:)); % azimuth and inclinations

    if isfield(protocol,'smalldel')
        M = size(protocol.smalldel,2);       
    elseif isfield(protocol,'smalldel1')
        M = size(protocol.smalldel1,2);   
    else
         M = size(protocol.G,1);       
    end
    E = zeros(M,Ndir);
  
    for i = 1:Ndir
        params = [x pi/2-incl(i) az(i)];
        E(:,i) = fun(params,protocol);
    end
 
        E = mean(E,2);
  
    
    % Compute the Jacobian matrix
if(nargout>1)
    J = zeros(length(E),1);
    dx = protocol.pert;
    if nargin < 3          
        xpert = x;
        protocol.diff=1;   
        xpert(1) = xpert(1)*(1+dx);    
        Epert = AstroSticks(xpert,protocol);
        dEtdD = (Epert - E)/(xpert(1)*dx);     
        
        J(:,1) = dEtdD;        
    else % x_deriv is a vector containing 1 and 0  which specifies the parameters that enter the jacobian      
            if x_deriv(1) ~= 0                 
                xpert = x;
                protocol.diff=1;   
                xpert(1) = xpert(1)*(1+dx);    
                Epert = AstroSticks(xpert,protocol,x_deriv);
                dEtdx = (Epert - E)/(xpert(1)*dx);
                J(:,1) = dEtdx;          
            end       
    end   
    protocol.diff=0;   
end   
    
end