function [ grp_delay , signal] = Onset_Det_n( FOUT,sampling_rate, smoothening_factor, winScaleFactor, thres )
    Mlf_file = 'result_1_ss.mlf';
    %FOUT=ones(length(FOUT),1).*abs(min(FOUT)) +FOUT +1;%*25;
    FOUT=abs(FOUT)+5;
   % winScaleFactor = 5:5:45;
    % Begin of paramterization
    winScaleFactor=10;
  
        for indexB = 1:length(smoothening_factor)
                for indexD = 1:length(thres)
                    %fprintf('Paramterization details : dwnsmpl_rate = %d\n smth_factor = %d, thres = %f\n',smoothening_factor(indexB),thres(indexD));
                    %======================================================================
                    % Part1: Using Amplitude Demodulation, and applying Group delay on it
                    S=FOUT;  
                    grp_delay = ones(length(S),1);
                    gd_sum = ones(length(S),1);

                    for wsfIndex = 1:length(winScaleFactor)
                        tempDir = sprintf('temp_%d',wsfIndex);
                        warning('off','all')
                        mkdir(tempDir); cd(tempDir);
                        energy_file_name = strcat(sprintf('neuron'),'.en');
                        dlmwrite(energy_file_name,S,'\n');
                        display(size(S));
                        spec_file_name = energy_file_name(1:end-2);
                        spec_file_name =strcat(spec_file_name,'spec');

                        % Invoking the binary                        
                        copyfile('../fe-words.base','fe-words.base');
                        ctrl_file = 'fe-words.base';
                        temp_ctrl_file = strcat('temp.base');

                        % Changing the winscalefactor parameter in config file
                        a = importdata(ctrl_file);
                        a = struct2cell(a);
                        a{1}(3) = winScaleFactor(wsfIndex);
                        %fprintf('Window scale factor is %d\n',winScaleFactor(wsfIndex));
                        fid0 = fopen(temp_ctrl_file,'w');
                        for i = 1:length(a{1})
                            fprintf(fid0,'%s %s %f\n',char(a{2}(i,1)),char(a{2}(i,2)),a{1}(i));
                        end
                        copyfile(temp_ctrl_file,ctrl_file);
                        %delete(temp_ctrl_file);
                        fclose(fid0);
                        %disp('here');    
                        dummy1 = 'b';
                        dummy2 = 'c';
                        dummy3 = 'd';
                        dummy4 = 'e';
                        dump = 'dump.txt';

                        system(sprintf('../WordSegmentWithSilenceRemoval %s %s %s %s %s %s %s > %s 2>&1',ctrl_file,energy_file_name,spec_file_name,dummy1,dummy2,dummy3,dummy4,dump));
                           %disp('here');
                        delete(energy_file_name);
                        temp = load(spec_file_name);
%			temp
%			display(size(grp_delay));
                        %delete(spec_file_name);
                        temp = temp(:,5);
                        temp(length(S)+1:end) = [];
                        %disp('here');
                        %figure;
                        %plot(temp);
                        grp_delay = grp_delay.*temp;
                        temp = temp - mean(temp);
                        gd_sum = gd_sum + cumsum(temp);
                        %disp('here');
                        cd ..; 
                    end

                    grp_delay = diff(gd_sum);
                    %grp_delay = temp;

                    %grp_delay = smooth(grp_delay,2*smoothening_factor(indexB),'moving');   % A moving average with 1 ms interval

                    grp_delay = grp_delay/max(grp_delay); 
                    grp_delay=[grp_delay(1:end-20);grp_delay(end-20).*ones(20,1)];
                    grp_delay=smooth(grp_delay,5);
                     %grp_delay = (grp_delay-mean(grp_delay))./std(grp_delay);
                    assignin('base','grp_delay',grp_delay);
                    %figure;
                    %======================================================================



                    %======================================================================
                    % Part2: Reading the contents of group delay file, and getting the
                    % onsets
                    threshold = thres(indexD);

                    stroke_loc = zeros(1,length(grp_delay));
                    % Go to each minima, and calculate height till next maxima. Keep a
                    % threshold on this to decide if stro
%	             display('GD dimension before stroke_loc');
%grp_delay;
                    t = 1:length(grp_delay);
                    [ymax,imax,ymin,imin] = extrema(grp_delay);

                    % sort the minimas and maximas;
 %                    imin
  %                   ymin
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


                     %==================================================================
                     % Algorithm1  for stroke location
                     index_stroke = 1;
                     peak_valley_heights = ymax - ymin;
                     peak_valley_heights = peak_valley_heights(1:length(peak_valley_heights));

                     for index = 1:1:length(peak_valley_heights)
                         if (peak_valley_heights(index) > threshold)
                             %stroke_loc(index_stroke) = ceil((imin(index) + imax(index))/2);
                             stroke_loc(index_stroke) = imin(index) ;
                             index_stroke = index_stroke + 2;
                         end
                     end
                     %==================================================================


                    stroke_loc(stroke_loc==0) = [];

                    assignin('base','stroke_loc',stroke_loc);
                    assignin('base','peaks',peaks);

                    %======================================================================



                    %======================================================================
                    % Printing in standard MLF format
                    dangerflag = 0;
                    X=S;
                    Fs=sampling_rate;

                    fid3 = fopen(Mlf_file,'a');

                    length_wav_file = length(X)*1/Fs;
                    stroke_loc = stroke_loc*sampling_rate/Fs; 
                    time=length(S)/Fs;% Converting into seconds
                    gdans=zeros(length(S),1);
                    gdans(stroke_loc)=1;
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
%		    display(size(signal));             
                end
        end
        
  
end

