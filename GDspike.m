%    Copyright (C) 2017 Speech and Music Technology Lab,
%     Indian Institute of Technology Madras
                
%     Contributed by Jilt Sebastian <jilt@cse.iitm.ac.in>
%    This file is a part of GDspike:Spike estimation evaluation system
%     GDspike is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or   
%    (at your option) any later version.

%    GDspike is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.                             

%    You should have received a copy of the GNU General Public License.          
%    If not, see <http://www.gnu.org/licenses/>.


%	 input parameters of the function GD_spike
%	 %% Load the input and pass it to the function: Ca_signal
%	 Threshold:          To be used in triangulation step (default= 0.4)
%	 Resampling_rate:    Rate of sampling used for obtaining the spike estimates (default=100)
%	 winScaleFactor:     GD_window_scale_factor- Division factor at the GD domain (default=4)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [ grp_delay, signal, GD_out] = GDspike(Ca_signal)

       
       %[Ca_signal]= load(Input); % signal and its sampling rate, edit appropriately
       
       %edit here
       Ca_Samp_rate    = 60;
       Resampling_rate = 60;
       threshold       = 0.4;
       winScaleFactor  = 4;

       %resampling the data
       resampled_Ca=resample(double(Ca_signal),Resampling_rate,Ca_Samp_rate);
       [Ca_signal_value]=resampled_Ca(1:floor(length(Ca_signal)/Ca_Samp_rate*Resampling_rate));
       [Ca_signal_time]= 0:1/Resampling_rate:length(Ca_signal_value)/Resampling_rate;       

 	S = Ca_signal_value;
	
        
	% GD computation======================================================================
	    grp_delay = ones(length(S),1);
        gd_sum = ones(length(S),1);

        tempDir = sprintf('temp_%d',winScaleFactor);
        warning('off','all')
        mkdir(tempDir); cd(tempDir);
        energy_file_name = strcat(sprintf('neuron'),'.en');
        dlmwrite(energy_file_name,S,'\n');
        %display(size(S));
        spec_file_name = energy_file_name(1:end-2);
        spec_file_name =strcat(spec_file_name,'spec');

        % Invoking the binary                        
        copyfile('../fe-words.base','fe-words.base');
        ctrl_file = 'fe-words.base';
        temp_ctrl_file = strcat('temp.base');

        % Changing the winscalefactor parameter in config file

        a = importdata(ctrl_file);
        a = struct2cell(a);
	    a{1}(3) = winScaleFactor;
        % fprintf('Window scale factor is %d\n',winScaleFactor(wsfIndex));
        fid0 = fopen(temp_ctrl_file,'w');
        for i = 1:length(a{1})
            fprintf(fid0,'%s %s %f\n',char(a{2}(i,1)),char(a{2}(i,2)),a{1}(i));
        end
        copyfile(temp_ctrl_file,ctrl_file);
        fclose(fid0);
   
        dummy1 = 'b';
        dummy2 = 'c';
        dummy3 = 'd';
        dummy4 = 'e';
        dump = 'dump.txt';
	
	% Part 1: Running the binary file-which gives the gd domain output

        system(sprintf('../WordSegmentWithSilenceRemoval %s %s %s %s %s %s %s > %s 2>&1',ctrl_file,energy_file_name,spec_file_name,dummy1,dummy2,dummy3,dummy4,dump));
        delete(energy_file_name);
        temp = load(spec_file_name);
        temp = temp(:,5);
        temp(length(S)+1:end) = [];
        grp_delay = grp_delay.*temp;
        temp = temp - mean(temp);
        gd_sum = gd_sum + cumsum(temp);
        cd ..; 
	    grp_delay = diff(gd_sum);
        grp_delay = grp_delay/max(grp_delay); 
        grp_delay=[grp_delay(1:end-20);grp_delay(end-20).*ones(20,1)];
        grp_delay=smooth(grp_delay,5);
        %grp_delay = (grp_delay-mean(grp_delay))./std(grp_delay);
        assignin('base','grp_delay',grp_delay);
      
        % Part2: Reading the contents of group delay file, and getting the spike locations
        
        spike_loc = zeros(1,length(grp_delay));
        t = 1:length(grp_delay);
        [ymax,imax,ymin,imin] = extrema(grp_delay);

        % sort the minimas and maximas;=================================================

        temp_min = sortrows([imin ymin]);
        imin = temp_min(:,1)';
        ymin = temp_min(:,2)';

        clear temp_min;

        temp_max = sortrows([imax ymax]);
        imax = temp_max(:,1)';
        ymax = temp_max(:,2)';
        clear temp_max;


      if (imin(1) < imax(1) )  % fine, just truncate the maximum
           imin(1) = []; ymin(1) = [];

 	      if (length(imin) > length(imax) )
 		   imin(length(imax)+1:end) = [];
           	   ymin(length(imax)+1:end) = [];
       	      elseif (length(imin) < length(imax) )
           	   imax(length(imin)+1:end) = [];
           	   ymax(length(imin)+1:end) = [];
       	      end
      else                                                    

           if (length(imin) > length(imax) )
               disp('this shouldnt have come');
               imin(length(imax)+1:end) = [];
               ymin(length(imax)+1:end) = [];
           elseif (length(imin) < length(imax) )
               imax(length(imin)+1:end) = [];
               ymax(length(imin)+1:end) = [];
           end
      end


      assignin('base','ymax',ymax);
      assignin('base','imax',imax);
      assignin('base','ymin',ymin);
      assignin('base','imin',imin);
      assignin('base','grp_delay',grp_delay);

     % A basic algorithm based on the threshold =============================================

     index_spike = 1;
     peak_valley_heights = ymax - ymin;
     peak_valley_heights = peak_valley_heights(1:length(peak_valley_heights));

     for index = 1:1:length(peak_valley_heights)
         if (peak_valley_heights(index) > threshold)
             spike_loc(index_spike) = imin(index) ;
             index_spike = index_spike + 2;
         end
     end
    %==================================================================

    % Part 3: Triangulation step

    spike_loc(spike_loc==0) = [];

    assignin('base','spike_loc',spike_loc);
    assignin('base','peaks',peaks);

    dangerflag = 0;
    X=S;
    Fs=Ca_Samp_rate;
    length_wav_file = length(X)*1/Fs;
    spike_loc = spike_loc*Ca_Samp_rate/Fs; 
    time=length(S)/Fs;% Converting into seconds
    gdans=zeros(length(S),1);
    gdans(spike_loc)=1;
    isort = sort([imax imin]);
    exts=0;
    exte=0;
    iexts=1;
    iexte=1;
    ind=1;

    for ind= 1:length(isort)-1
        iexts=isort(ind);
        iexte=isort(ind+1);
        exts = grp_delay(iexts);
        exte = grp_delay(iexte);
        valley_lenght=abs(exts-exte);
        signal(iexts)=0;
        signal(iexte)=0;
        intermediate_samples=iexte-iexts-2;
        increment_value=valley_lenght/(floor((intermediate_samples+1)/2));
        sum=0;
        if valley_lenght>0.1*max(grp_delay)
            for j=iexts+1:iexts+floor((intermediate_samples+1)/2)
                sum=sum+increment_value;
                signal(j)=sum;
            end 
            sum=0;
            for j=iexte-1:-1:iexts+floor((intermediate_samples+1)/2)+1
                sum=sum+increment_value;
                signal(j)=sum;
            end
        else
            signal(iexts:iexte)=0;
        end
    end
    signal(isort(length(isort)):time*Fs)=0;


    % end of GD computation====================================================================== 
    
   GD_out=downsample_(Inc_Spk1(signal',threshold),1);
end


