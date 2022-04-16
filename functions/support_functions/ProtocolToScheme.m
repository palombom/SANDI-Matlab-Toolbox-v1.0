function ProtocolToScheme(protocol, schemefile,tau)
%
% camino.m--------------------------------------------------------------
% Write a schemefile from a protocol
% 
% ProtocolToScheme(protocol, schemefile)
% 
% Description: Writes a Camino Version 1 or GRADIENT_WAVEFORM schemefile 
% from protocol object
%
% Parameters:  
% protocol - structure containing diffusion acquision information 
% schemefile -  Camino Version 1/GRADIENT_WAVEFORM schemefile
%   
%------------------------------------------------------------------------
% This file is part of the camino.m toolbox.
% Copyright (c) 2015, UCL Microstructure Imaging Group (MIG), All rights reserved.
% Distributed under the Modified BSD Licence (see: LICENSE.pdf).
%
% Authors:
%   Andrada Ianus (a.ianus.11@ucl.ac.uk)
%   Daniel C Alexander (d.alexander@ucl.ac.uk)
%   Gary Hui Zhang (gary.zhang@ucl.ac.uk)
%  

fid = fopen(schemefile, 'w', 'b');

if(strcmp(protocol.pulseseq, 'PGSE')) || (strcmp(protocol.pulseseq, 'STEAM'))
    % Write in the header
    fprintf(fid, 'VERSION: 1\n');
   display(['Writing schemefile for ' protocol.pulseseq]);
   % Write out each measurement
   if isfield(protocol,'totalmeas')
        M = protocol.totalmeas;
    else
        M = length(protocol.smalldel);
    end
    if(strcmp(protocol.pulseseq, 'PGSE'))
      TE = max(protocol.delta+protocol.smalldel);
    elseif(strcmp(protocol.pulseseq, 'STEAM'))
      TE = max(2*protocol.smalldel);
    end
    for i=1:M
      if(strcmp(protocol.testrategy, 'variable'))
        if(strcmp(protocol.pulseseq, 'PGSE'))
          TE=protocol.delta(i)+protocol.smalldel(i);
        elseif(strcmp(protocol.pulseseq, 'STEAM'))
          TE = 2*protocol.smalldel(i);
        end
      elseif (strcmp(protocol.testrategy, 'given'))
          
          TE = protocol.TE(i);
      end
      fprintf(fid, '%f\t%f\t%f\t%f\t%f\t%f\t%f\t\n', protocol.grad_dirs(i,1), protocol.grad_dirs(i,2), protocol.grad_dirs(i,3), protocol.G(i), protocol.delta(i), protocol.smalldel(i), TE);
    end
else
    % Write in the header
    fprintf(fid, 'VERSION: GRADIENT_WAVEFORM\n');

    if isfield(protocol,'tau')
        kdel = protocol.tau;
    else
        kdel = tau;
    end
    if isfield(protocol,'totalmeas')
        M = protocol.totalmeas;
    else
        M = length(protocol.smalldel);
    end
    if(strcmp(protocol.pulseseq, 'FullSTEAM'))
        display(['Writing schemefile for ' protocol.pulseseq]);
        ks = 0:kdel:(max(protocol.TE+protocol.TM)+(kdel*100));
        for j=1:M
            wav = zeros(3,length(ks));     
            ddel = protocol.smalldel(j);
            dG = protocol.G(j).*protocol.grad_dirs(j,:);
        
            for i=1:length(ks)
                if(ks(i)<protocol.TE(j)/2-protocol.smalldel(j)-protocol.gap1(j)-protocol.sdelc(j)-protocol.sdelr(j))
                    wav(:,i) = 0;
                elseif(ks(i)<protocol.TE(j)/2-protocol.gap1(j)-protocol.sdelc(j)-protocol.sdelr(j))
                    wav(:,i) = dG;
                elseif(ks(i)<protocol.TE(j)/2-protocol.sdelc(j)-protocol.sdelr(j))
                    wav(:,i) = 0;
                elseif(ks(i)<protocol.TE(j)/2-protocol.sdelr(j))
                    wav(:,i) = protocol.cG(j,:);
                elseif(ks(i)<protocol.TE(j)/2)
                    wav(:,i) = protocol.rG(j,:);
                elseif(ks(i)<protocol.TE(j)/2+protocol.TM(j))
                    wav(:,i) = 0;
                elseif(ks(i)<protocol.TE(j)/2+protocol.sdelr(j)+protocol.TM(j))
                    wav(:,i) = -protocol.rG(j,:);
                elseif(ks(i)<protocol.TE(j)/2+protocol.sdelc(j)+protocol.sdelr(j)+protocol.TM(j))
                    wav(:,i) = -protocol.cG(j,:);
                elseif(ks(i)<protocol.TE(j)/2+protocol.sdelc(j)+protocol.sdelr(j)+protocol.gap2(j)+protocol.TM(j))
                    wav(:,i) = 0;
                elseif(ks(i)<protocol.TE(j)/2+protocol.sdelc(j)+protocol.sdelr(j)+protocol.gap2(j)+protocol.smalldel(j)+protocol.TM(j))
                    wav(:,i) = -dG;
                elseif(ks(i)<protocol.TE(j)+protocol.TM(j))
                    wav(:,i) = 0;
                end
            end  
            % Construct the output string
        str = sprintf('%i %f', length(ks), kdel);
        for i=1:length(ks)
            str = sprintf('%s %f %f %f', str, wav(1,i), wav(2,i), wav(3,i));
        end
        fprintf(fid, '%s\n', str);
        end    
        
    elseif (strcmp(protocol.pulseseq, 'OGSE'))
        display(['Writing schemefile for ' protocol.pulseseq]);
      if(isfield(protocol,'omega'))
         total_time = max(ceil((protocol.smalldel+protocol.delta)/kdel)).*kdel+kdel;
         ks = 0:kdel:total_time;
         wav = zeros(3,length(ks));
        for j = 1:M
        dG =protocol.G(j).*protocol.grad_dirs(j,:);
            if ~isfield(protocol,'phase')  % no phase            
                for i=2:length(ks)
                    if(ks(i)-protocol.smalldel(j)<-1E-10)       
                        wav(:,i) =dG*sin(protocol.omega(j)*ks(i));
                        i_smalldel = i+1;
                        i_delta = i_smalldel; % for the case when delta-smalldel = 0; needed  for some tests
                    elseif(ks(i)-protocol.delta(j)<-1E-10)    
                        wav(:,i) = 0;
                        i_delta = i;
                   elseif i-i_delta <=i_smalldel % operate on integers
                         if ~isfield(protocol,'mirror')  || protocol.mirror ==0 % repeated
                            wav(:,i) = -wav(:,i-i_delta);
                         else
                            wav(:,i) = -wav(:,i_smalldel-(i-i_delta)+1);
                         end                      
                    else
                        wav(:,i) = 0;
                    end  
                end            
            else           
                for i=2:length(ks)                     
                    if(ks(i)-protocol.smalldel(j)<-1E-10)
                          wav(:,i) =  dG*sin(protocol.omega(j)*ks(i)-protocol.phase(j));
                          i_smalldel = i+1;
                          i_delta = i_smalldel; % for the case when delta-smalldel = 0; needed  for some tests
                    elseif(ks(i)-protocol.delta(j)<-1E-10)
                          wav(:,i) = 0;
                          i_delta =  i;
                    elseif i-i_delta <=i_smalldel % operate on integers                       
                         if  ~isfield(protocol,'mirror') || protocol.mirror ==0 % repeated
                            wav(:,i) = -wav(:,i-i_delta);
                         else
                            wav(:,i) = -wav(:,i_smalldel-(i-i_delta)+1);
                         end                       
                    else
                         wav(:,i) = 0;
                    end
                end           
            end        
            % Construct the output string
            str = sprintf('%i %f', length(ks), kdel);
            for i=1:length(ks)
                str = sprintf('%s %f %f %f', str, wav(1,i), wav(2,i), wav(3,i));
            end
            fprintf(fid, '%s\n', str);
        end    
     else
        error('the protocol does not contain enough information to generate SWOGSE');
     end
        
    elseif (strcmp(protocol.pulseseq, 'SWOGSE'))
         display(['Writing schemefile for ' protocol.pulseseq]);
     if(isfield(protocol,'omega'))
         total_time = max(ceil((protocol.smalldel+protocol.delta)/kdel)).*kdel+kdel;
         ks = 0:kdel:total_time;
         wav = zeros(3,length(ks));
        for j = 1:M
        dG =protocol.G(j).*protocol.grad_dirs(j,:);
            if ~isfield(protocol,'phase')  % no phase            
                for i=2:length(ks)
                    if(ks(i)-protocol.smalldel(j)<-1E-10)       
                        wav(:,i) =dG*(-1).^floor(protocol.omega(j)./pi()*ks(i)-1E-10);
                        i_smalldel = i+1;
                        i_delta = i_smalldel; % for the case when delta-smalldel = 0; needed  for some tests
                    elseif(ks(i)-protocol.delta(j)<-1E-10)    
                        wav(:,i) = 0;
                        i_delta = i;
                   elseif i-i_delta <=i_smalldel % operate on integers
                         if ~isfield(protocol,'mirror')  || protocol.mirror ==0 % repeated
                            wav(:,i) = -wav(:,i-i_delta);
                         else
                            wav(:,i) = -wav(:,i_smalldel-(i-i_delta)+1);
                         end                      
                    else
                        wav(:,i) = 0;
                    end  
                end            
            else           
                for i=2:length(ks)                     
                    if(ks(i)-protocol.smalldel(j)<-1E-10)
                          wav(:,i) =  dG*(-1).^floor((protocol.omega(j)*ks(i)-protocol.phase(j))./pi()-1E-10);
                          i_smalldel = i+1;
                          i_delta = i_smalldel; % for the case when delta-smalldel = 0; needed  for some tests
                    elseif(ks(i)-protocol.delta(j)<-1E-10)
                          wav(:,i) = 0;
                          i_delta =  i;
                    elseif i-i_delta <=i_smalldel % operate on integers                       
                         if  ~isfield(protocol,'mirror') || protocol.mirror ==0 % repeated
                            wav(:,i) = -wav(:,i-i_delta);
                         else
                            wav(:,i) = -wav(:,i_smalldel-(i-i_delta)+1);
                         end                       
                    else
                         wav(:,i) = 0;
                    end
                end           
            end        
            % Construct the output string
            str = sprintf('%i %f', length(ks), kdel);
            for i=1:length(ks)
                str = sprintf('%s %f %f %f', str, wav(1,i), wav(2,i), wav(3,i));
            end
            fprintf(fid, '%s\n', str);
        end    
     else
        error('the protocol does not contain enough information to generate SWOGSE');
     end     
    elseif (strcmp(protocol.pulseseq, 'TWOGSE'))
         display(['Writing schemefile for ' protocol.pulseseq]);
     if(isfield(protocol,'slew_rate') && isfield(protocol,'omega') )
        total_time = max(ceil((protocol.smalldel+protocol.delta)/kdel)).*kdel+kdel;
        ks = 0:kdel:total_time;
        wav = zeros(3,length(ks));  
        for j = 1:M            
            rt = protocol.G(j)/protocol.slew_rate(j); % rise time  
            Nt = floor(protocol.smalldel(j).*protocol.omega(j)./pi+0.00000000001);
            dG =protocol.G(j).*protocol.grad_dirs(j,:);
            if (protocol.smalldel(j)-Nt.*pi/protocol.omega(j)>2*rt) 
                for i=2:length(ks)            
                    if(ks(i)<Nt.*pi/protocol.omega(j))
                        it = floor(protocol.omega(j)*ks(i)./pi()+0.00000001);
                        if( ks(i)<it *pi/protocol.omega(j)+rt)
                            wav(:,i) = dG./rt.*(-1).^it.*(ks(i)-it.*pi./protocol.omega(j));
             
                        elseif( ks(i)<(it+1)*pi/protocol.omega(j)-rt)
                            wav(:,i) = dG*(-1)^it;
                        else                   
                            wav(:,i) = dG./rt.*(-1).^it.*((it+1).*pi./protocol.omega(j)-ks(i));
                        end
                    elseif(ks(i)<protocol.smalldel(j)) % last oscillation
                        if( ks(i)<Nt.*pi/protocol.omega(j)+rt)
                            wav(:,i) = dG./rt.*(-1).^Nt.*(ks(i)-Nt.*pi./protocol.omega(j));
                 
                        elseif( ks(i)<protocol.smalldel(j)-rt)
                            wav(:,i) = dG*(-1)^Nt;
                        else
                            wav(:,i) = dG./rt.*(-1).^Nt.*(protocol.smalldel(j)-ks(i));
                        end
                        i_smalldel = i;
                    elseif(ks(i)<protocol.delta(j))
                        wav(:,i) = 0;
                        i_delta = i;
                    elseif i-i_delta <=i_smalldel % operate on integers                        
                         if ~isfield(protocol,'mirror')  || protocol.mirror ==0 % repeated
                            wav(:,i) = -wav(:,i-i_delta);
                         else
                            wav(:,i) = -wav(:,i_smalldel-(i-i_delta)+1);
                         end
                    else
                        wav(:,i) = 0;
                    end
                end
            else % case of very small last oscillation or integer half periods
                for i=2:length(ks)            
                    if(ks(i)<Nt.*pi/protocol.omega(j))
                        it = floor(protocol.omega(j)*ks(i)./pi()+0.00000001);
                        if( ks(i)<it *pi/protocol.omega(j)+rt)
                            wav(:,i) = dG./rt.*(-1).^it.*(ks(i)-it.*pi./protocol.omega(j));
                        elseif( ks(i)<(it+1)*pi/protocol.omega(j)-rt)
                            wav(:,i) = dG*(-1)^it;
                        else                  
                            wav(:,i) = dG./rt.*(-1).^it.*((it+1).*pi./protocol.omega(j)-ks(i));                    
                        end
                        i_smalldel = i;
                    elseif(ks(i)<protocol.smalldel(j)) % last oscillation
                        if( ks(i)<Nt.*pi/protocol.omega(j)+(protocol.smalldel(j)-Nt.*pi/protocol.omega(j))/2)
                            wav(:,i) = dG./rt.*(-1).^Nt.*(ks(i)-Nt.*pi./protocol.omega(j));  
                        else
                            wav(:,i) = dG./rt.*(-1).^Nt.*(protocol.smalldel(j)-ks(i));
                        end
                        i_smalldel = i;
                    elseif(ks(i)<protocol.delta(j))
                        wav(:,i) = 0;
                        i_delta = i;
                    elseif i-i_delta <=i_smalldel % operate on integers                       
                         if ~isfield(protocol,'mirror')  || protocol.mirror ==0 % repeated
                            wav(:,i) = -wav(:,i-i_delta);
                         else
                            wav(:,i) = -wav(:,i_smalldel-(i-i_delta)+1);
                         end                            
                    else
                        wav(:,i) = 0;
                    end
                end
            end
            % Construct the output string
            str = sprintf('%i %f', length(ks), kdel);
            for i=1:length(ks)
                str = sprintf('%s %f %f %f', str, wav(1,i), wav(2,i), wav(3,i));
            end
            fprintf(fid, '%s\n', str);
        end
     else
        error('the protocol does not contain enough information to generate TWOGSE');
     end   
    elseif (strcmp(protocol.pulseseq, 'GEN'))
        display(['Writing schemefile for ' protocol.pulseseq]);
        Gx = protocol.G(:,1:3:end);
        Gy = protocol.G(:,2:3:end);
        Gz = protocol.G(:,3:3:end); 
        for j = 1:size(protocol.G,1)
            % Construct the output string
            str = sprintf('%i %f', size(Gx,2), kdel);
            for i=1:size(Gx,2)
                str = sprintf('%s %f %f %f', str, Gx(j,i), Gy(j,i), Gz(j,i));
            end
            fprintf(fid, '%s\n', str);
        end
    else  error(['Not implemented for pulse sequence: ', protocol.pulseseq]);
    end        

end
  
fclose(fid);