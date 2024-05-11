%test = edfread('C:\Users\Pragati\Desktop\PhD_Cognitive Science\Semester_4\Neuroinformatics_spring_2024\Course_Project\Acute_stroke_\edffile\edffile\sub-01\eeg\sub-01_task-motor-imagery_eeg.edf');
%t = readtable('C:\Users\Pragati\Desktop\PhD_Cognitive Science\Semester_4\Neuroinformatics_spring_2024\Course_Project\Acute_stroke_\task-motor-imagery_events.tsv', "FileType","text",'Delimiter', '\t');
% Pre-processing for all patients (left and right hemiplegia)
participant_id_left = [2,5,6,8,10,18,21,40,49,50];
participant_id_right = [1,7,12,14, 32,33,38,41,45];

participant_id = participant_id_right;

for l = 1:1:length(participant_id)
    k = participant_id(l);
    if k <10
        load(sprintf('C:\\Users\\Pragati\\Desktop\\PhD_Cognitive Science\\Semester_4\\Neuroinformatics_spring_2024\\Course_Project\\Acute_stroke_\\sourcedata\\sourcedata\\sub-0%d\\sub-0%d_task-motor-imagery_eeg.mat', k, k));
    else
        load(sprintf('C:\\Users\\Pragati\\Desktop\\PhD_Cognitive Science\\Semester_4\\Neuroinformatics_spring_2024\\Course_Project\\Acute_stroke_\\sourcedata\\sourcedata\\sub-%d\\sub-%d_task-motor-imagery_eeg.mat', k, k));
    end
  permuted_eeg = permute(eeg.rawdata, [2 3 1]); % permuting into new dims = channels * timepoints * trials 
    % We recognized 33 as marker channel and 18 as ref channel - dropping them
    % alongside last 2 EOG channels
    permuted_eeg_drop = [permuted_eeg(1:17,:,:); permuted_eeg(19:30,:,:)];
    %trial_data = squeeze(permuted_eeg_drop(:,:,1)); %We verified the drop
    eeglab;
    EEG = pop_importdata('dataformat','array','nbchan',size(permuted_eeg_drop,1),'data','permuted_eeg_drop','setname','sub_0','srate',500,'pnts',size(permuted_eeg_drop,2),'chanlocs','C:\\Users\\Pragati\\Desktop\\PhD_Cognitive Science\\Semester_4\\Neuroinformatics_spring_2024\\Course_Project\\Acute_stroke_\\task-motor-imagery_electrodes - Copy.tsv');
    EEG.times = linspace(0,8,4000);
    EEG.data = EEG.data - repmat(mean(EEG.data,2),[1 EEG.pnts 1]); % Remove DC offset
% % time frequency decomposition plot for all data channels (before preprocessing) 
%     figure; 
%     metaplottopo( EEG.data, 'plotfunc', 'newtimef', 'chanlocs', EEG.chanlocs, 'plotargs', ...
%                    {EEG.pnts, [EEG.xmin EEG.xmax]*1000, EEG.srate, [0], 'plotitc', 'off', 'ntimesout', 50, 'padratio', 1});
% spectopo plot (before preprocessing) 
    figure;
    spectopo(EEG.data,EEG.pnts,EEG.srate); 
    % Butterfly Plot (before pre-processing)
     x_axis_limit = size(EEG.data,2); % in ms
     figure;
     for chan= 1:26
         plot(1: x_axis_limit,squeeze(mean(EEG.data (chan,:,:),3)))
         title(["Butterfly plot pre processing for partcipant" num2str(k)])
         hold on; 
      end
    % EEG.data = detrend(EEG.data);
    %EEG = pop_eegfiltnew(EEG, 'locutoff',0.1); % High passfilter
    % EEG = pop_eegfiltnew(EEG, 'locutoff', 0.5, 'hicutoff', 51);
    EEG = pop_eegfiltnew(EEG, 'locutoff', 49, 'hicutoff', 51, 'revfilt', 1);%band-stop filtering
    EEG = pop_eegfiltnew(EEG, 1, 0); %High pass filter

    eeglab redraw;
    %Average Rereferncing 
    for i = 1:EEG.trials    
         y = EEG.data(:,:,i);
         reference_y = mean(y,1);
         reference_y_stacked = repmat(reference_y,29,1);
         reference_y_referenced(:,:,i) = y - reference_y_stacked;
    end
    EEG.data = reference_y_referenced;
    figure;
    spectopo(EEG.data,EEG.pnts,EEG.srate); 
    % Artifact removal
    permuted_eeg_after_rereference = permute(EEG.data, [1,3,2]);
    [iChanKeep,iEvKeep,~] = jwcleaneegevents_v01(permuted_eeg_after_rereference, 1, 'artifactremoval', [2.3 2.3]);
    permuted_eeg_after_rereference_drop = permuted_eeg_after_rereference(iChanKeep,iEvKeep,:);
    permuted_eeg_after_rereference_final = permute(permuted_eeg_after_rereference_drop, [1,3,2]);
    EEG.data = permuted_eeg_after_rereference_final;
    eeglab redraw;
    
%     clean_rawdata(EEG)- for ICA
%     [~, x] = clean_channels(EEG);
%     [ALLEEG, EEG,CURRERNTSET] = processMARA(ALLEEG,EEG,CURRENTSET,[0,0,0,0,0]) ;% Performing Runica before MARA classification 
    % Select and remove MARA rejected components
%     EEG = pop_selectcomps_MARA(EEG, EEG.reject.gcompreject); %select components for rejection
%     rejected_indices = find (EEG.reject.gcompreject);
%     EEG = pop_subcomp(EEG,rejected_indices,0); %rejects components from indices marked for rejection

% % Time-frequnency decomposition plot for all data channels (after preprocessing)  
%     figure; 
%     metaplottopo( EEG.data, 'plotfunc', 'newtimef', 'chanlocs', EEG.chanlocs, 'plotargs', ...
%                    {EEG.pnts, [EEG.xmin EEG.xmax]*1000, EEG.srate, [0], 'plotitc', 'off', 'ntimesout', 50, 'padratio', 1});
%    
% spectopo plot (post- preprocessing) 
 figure;
    spectopo(EEG.data,EEG.pnts,EEG.srate);
% Butterfly Plots (post-processing) 
figure;
     for chan= 1:size(EEG.data,1)
         plot(1: x_axis_limit,squeeze(mean(EEG.data (chan,:,:),3)))
         title(["Butterfly plot post processing for partcipant" num2str(k)])
         hold on; 
      end
% saving the output 
    save(fullfile("C:\Users\Pragati\git\Neuroinformatics\Course_Project_Stroke\Right_paralysis_preprocessed", "EEG_processed" + num2str(k) + ".mat"),"EEG");
    save(fullfile("C:\Users\Pragati\git\Neuroinformatics\Course_Project_Stroke\Right_paralysis_preprocessed", "EEG_processed" + num2str(k) + ".set"),"EEG");
    save (fullfile("C:\Users\Pragati\git\Neuroinformatics\Course_Project_Stroke\Right_paralysis_preprocessed", "EEG_dynamic_trials" + num2str(k) + ".mat"),"iEvKeep")
%  close all;
    eeglab redraw; 

    close all;
end

%% Extracting epochs (3 stages- Instruction, Motor Imagery, Break)

for l = 1:1:length(participant_id)
    k = participant_id(l);
    dynamic_trial= load (['C:\Users\Pragati\git\Neuroinformatics\Course_Project_Stroke\Right_paralysis_preprocessed\' 'EEG_dynamic_trials' num2str(k) '.mat']);
    eve = dynamic_trial.iEvKeep;
    eeglab;
    close all;
    filename = ['EEG_processed' num2str(k) '.set'];
    EEG = pop_loadset('filename',filename,'filepath','C:\\Users\\Pragati\\git\\Neuroinformatics\\Course_Project_Stroke\\Right_paralysis_preprocessed');
    epoch_lengths = [2, 4, 2];% lengths of each epochs of interest in seconds
    epoch_lengths_points = epoch_lengths * EEG.srate;
    num_channels = EEG.nbchan;
    num_trials = EEG.trials;

    inst_epoch = zeros(num_channels, epoch_lengths_points(1), num_trials);
    imagery_epoch_l1 = zeros(num_channels, epoch_lengths_points(2), num_trials);
    break_epoch_l1 = zeros(num_channels, epoch_lengths_points(3), num_trials);
    
    l1 = 1; %for left hand motor imagery 
    l2 = 1; %for right hand motor imagery

    for m = 1:EEG.trials  
        trial_data = EEG.data (:, :, m);
        inst_epoch(:, :, m) = trial_data(:, 1:epoch_lengths_points(1));
        imagery_epoch_l1(:, :, m) = trial_data(:, (epoch_lengths_points(1)+1):(epoch_lengths_points(1)+epoch_lengths_points(2)));
        break_epoch_l1(:, :, m) = trial_data(:, (end-epoch_lengths_points(3)+1):end);
        event = mod(eve(m),2); 

        if event == 0
            l2_inst_epoch(:, :, l2) = trial_data(:, 1:epoch_lengths_points(1));
            l2_imagery_epoch(:, :, l2) = trial_data(:, (epoch_lengths_points(1)+1):(epoch_lengths_points(1)+epoch_lengths_points(2)));
            l2_break_epoch(:, :, l2) = trial_data(:, (end-epoch_lengths_points(3)+1):end);
            l2 = l2 + 1;
        else
            l1_inst_epoch(:, :, l1) = trial_data(:, 1:epoch_lengths_points(1));
            l1_imagery_epoch(:, :, l1) = trial_data(:, (epoch_lengths_points(1)+1):(epoch_lengths_points(1)+epoch_lengths_points(2)));
            l1_break_epoch(:, :, l1) = trial_data(:, (end-epoch_lengths_points(3)+1):end);
            l1 = l1 + 1;
        end
    end
% saving the output 
    save(fullfile("C:\Users\Pragati\git\Neuroinformatics\Course_Project_Stroke\Right_paralysis_preprocessed", "EEG_processed_epoched_l1_inst_epoch" + num2str(k) + ".mat"),"l1_inst_epoch");
    save(fullfile("C:\Users\Pragati\git\Neuroinformatics\Course_Project_Stroke\Right_paralysis_preprocessed", "EEG_processed_epoched_l1_inst_epoch" + num2str(k) + ".set"),"l1_inst_epoch");
    save(fullfile("C:\Users\Pragati\git\Neuroinformatics\Course_Project_Stroke\Right_paralysis_preprocessed", "EEG_processed_epoched_l1_imagery_epoch" + num2str(k) + ".set"),"l1_imagery_epoch");
    save(fullfile("C:\Users\Pragati\git\Neuroinformatics\Course_Project_Stroke\Right_paralysis_preprocessed", "EEG_processed_epoched_l1_imagery_epoch" + num2str(k) + ".mat"),"l1_imagery_epoch");
    save(fullfile("C:\Users\Pragati\git\Neuroinformatics\Course_Project_Stroke\Right_paralysis_preprocessed", "EEG_processed_epoched_l1_break_epoch" + num2str(k) + ".set"),"l1_break_epoch");
    save(fullfile("C:\Users\Pragati\git\Neuroinformatics\Course_Project_Stroke\Right_paralysis_preprocessed", "EEG_processed_epoched_l1_break_epoch" + num2str(k) + ".mat"),"l1_break_epoch");
    save(fullfile("C:\Users\Pragati\git\Neuroinformatics\Course_Project_Stroke\Right_paralysis_preprocessed", "EEG_processed_epoched_l2_inst_epoch" + num2str(k) + ".mat"),"l2_inst_epoch");
    save(fullfile("C:\Users\Pragati\git\Neuroinformatics\Course_Project_Stroke\Right_paralysis_preprocessed", "EEG_processed_epoched_l2_inst_epoch" + num2str(k) + ".set"),"l2_inst_epoch");
    save(fullfile("C:\Users\Pragati\git\Neuroinformatics\Course_Project_Stroke\Right_paralysis_preprocessed", "EEG_processed_epoched_l2_imagery_epoch" + num2str(k) + ".set"),"l2_imagery_epoch");
    save(fullfile("C:\Users\Pragati\git\Neuroinformatics\Course_Project_Stroke\Right_paralysis_preprocessed", "EEG_processed_epoched_l2_imagery_epoch" + num2str(k) + ".mat"),"l2_imagery_epoch");
    save(fullfile("C:\Users\Pragati\git\Neuroinformatics\Course_Project_Stroke\Right_paralysis_preprocessed", "EEG_processed_epoched_l2_break_epoch" + num2str(k) + ".set"),"l2_break_epoch");
    save(fullfile("C:\Users\Pragati\git\Neuroinformatics\Course_Project_Stroke\Right_paralysis_preprocessed", "EEG_processed_epoched_l2_break_epoch" + num2str(k) + ".mat"),"l2_break_epoch"); 

%     figure; % plotting time domain signal for one channel, one trial and one of the labels
%     subplot(3,1,1);
%     plot(inst_epoch(chan2plot,:,trial2plot));
%     title(['Time Domain "Instruction" signal for channel no. ' , num2str(chan2plot) , "Trial no. ",num2str(trial2plot) , " Event no. " , event]);
%     subplot(3,1,2);
%     plot(imagery_epoch(chan2plot,:,trial2plot));
%     title(['Time Domain "Motor Imagery" signal for channel no. ' , num2str(chan2plot) , "Trial no. ",num2str(trial2plot) , " Event no. " , event] );
%     subplot(3,1,3);
%     plot(break_epoch(chan2plot,:,trial2plot));
%     title(['Time Domain "Break" signal for channel no. ' , num2str(chan2plot) , "Trial no. ",num2str(trial2plot) , " Event no. " , event]);
% 
%     figure; % plotting frequency domain signal for one channel, one trial and one of the labels (fft)
%     subplot(3,1,1);
%     plot(abs(fft(inst_epoch(chan2plot,:,trial2plot))));
%     title(['Frequency Domain "Instruction" signal for channel no. ' , num2str(chan2plot) , "Trial no. ",num2str(trial2plot) , " Event no. " , event]);
%     subplot(3,1,2);
%     plot(abs(fft(imagery_epoch(chan2plot,:,trial2plot))));
%     title(['Frequency Domain "Motor Imagery" signal for channel no. ' , num2str(chan2plot) , "Trial no. ",num2str(trial2plot) , " Event no. " , event]);
%     subplot(3,1,3);
%     plot(abs(fft(break_epoch(chan2plot,:,trial2plot))));
%     title(['Frequency Domain "Break" signal for channel no. ' , num2str(chan2plot) , "Trial no. ",num2str(trial2plot) , " Event no. " , event]);
%     
  
  
    clear l1_inst_epoch l1_imagery_epoch l1_break_epoch l2_inst_epoch l2_imagery_epoch l2_break_epoch dynamic_trial;
%close all;
end


%% Averaging across all participants all trials for label 1 and label 2 
% We need to change the chan variable ie C3 and C4, and the change the
% location for the left and right hemiplegia respectively. 
l1_imag_avg = zeros(1,2000);
l2_imag_avg = zeros(1,2000);

chan = 14;
for l=1:1:length(participant_id)
    k = participant_id(l);
    instruction_epoch_l1 = load (['C:\Users\Pragati\git\Neuroinformatics\Course_Project_Stroke\Left_paralysis_preprocessed\' 'EEG_processed_epoched_l1_inst_epoch' num2str(k) '.mat']);
    l1_inst_epoch = instruction_epoch_l1.l1_inst_epoch;
    imagery_epoch_l1= load (['C:\Users\Pragati\git\Neuroinformatics\Course_Project_Stroke\Left_paralysis_preprocessed\' 'EEG_processed_epoched_l1_imagery_epoch' num2str(k) '.mat']);
    l1_imagery_epoch = imagery_epoch_l1.l1_imagery_epoch;
    break_epoch_l1 = load (['C:\Users\Pragati\git\Neuroinformatics\Course_Project_Stroke\Left_paralysis_preprocessed\' 'EEG_processed_epoched_l1_break_epoch' num2str(k) '.mat']);
    l1_break_epoch = break_epoch_l1.l1_break_epoch;
    instruction_epoch_l2 = load (['C:\Users\Pragati\git\Neuroinformatics\Course_Project_Stroke\Left_paralysis_preprocessed\' 'EEG_processed_epoched_l2_inst_epoch' num2str(k) '.mat']);
    l2_inst_epoch = instruction_epoch_l2.l2_inst_epoch;
    imagery_epoch_l2 = load (['C:\Users\Pragati\git\Neuroinformatics\Course_Project_Stroke\Left_paralysis_preprocessed\' 'EEG_processed_epoched_l2_imagery_epoch' num2str(k) '.mat']);
    l2_imagery_epoch = imagery_epoch_l2.l2_imagery_epoch;
    break_epoch_l2 = load (['C:\Users\Pragati\git\Neuroinformatics\Course_Project_Stroke\Left_paralysis_preprocessed\' 'EEG_processed_epoched_l2_break_epoch' num2str(k) '.mat']);
    l2_break_epoch = break_epoch_l2.l2_break_epoch;
   
%     % Plots of spectotopo for each of the 3 stages 
%     figure;
%     subplot(2,1,1);
%     spectopo(l1_inst_epoch,1000,EEG.srate); 
% 
%     subplot(2,1,2);
%     spectopo(l2_inst_epoch,1000,EEG.srate); 
% 
%     figure;
%     subplot(2,1,1);
%     spectopo(l1_imagery_epoch,2000,EEG.srate); 
% 
%     subplot(2,1,2);
%     spectopo(l2_imagery_epoch,2000,EEG.srate); 
% 
%     figure;
%     subplot(2,1,1);
%     spectopo(l1_break_epoch,1000,EEG.srate); 
% 
%     subplot(2,1,2);
%     spectopo(l2_break_epoch,1000,EEG.srate); 
    currsignal =  l1_imagery_epoch;

    channel_index = chan;
    x_axis_limit = size(currsignal,2); % in ms

    l1_imag_avg = l1_imag_avg + squeeze(mean(currsignal(channel_index,:,:),3)) / length(participant_id);

%     figure;
%     plot(1: x_axis_limit,squeeze(mean(currsignal(channel_index,:,:),3)))

    currsignal =  l2_imagery_epoch;
    channel_index = chan;
    x_axis_limit = size(currsignal,2); % in ms

    l2_imag_avg = l2_imag_avg + squeeze(mean(currsignal(channel_index,:,:),3)) / length(participant_id);

%      figure;
%      plot(1: x_axis_limit,squeeze(mean(currsignal(channel_index,:,:),3)))
end

figure;
subplot(211)
plot(l1_imag_avg);
title(" Averaged Imagery epoch across all Left Hemiplegia patients for label 1 for channel C3");
subplot(212)
plot(l2_imag_avg);
title("Averaged Imagery epoch across all Left Hemiplegia patients for label 2 for channel C3");


%% Wavelet Convolution

wav_amp = 1;

for l=1:1:1
    k = participant_id(l);
    instruction_epoch_l1 = load (['C:\Users\Pragati\git\Neuroinformatics\Course_Project_Stroke\Right_paralysis_preprocessed\' 'EEG_processed_epoched_l1_inst_epoch' num2str(k) '.mat']);
    l1_inst_epoch = instruction_epoch_l1.l1_inst_epoch;
    imagery_epoch_l1= load (['C:\Users\Pragati\git\Neuroinformatics\Course_Project_Stroke\Right_paralysis_preprocessed\' 'EEG_processed_epoched_l1_imagery_epoch' num2str(k) '.mat']);
    l1_imagery_epoch = imagery_epoch_l1.l1_imagery_epoch;
    break_epoch_l1 = load (['C:\Users\Pragati\git\Neuroinformatics\Course_Project_Stroke\Right_paralysis_preprocessed\' 'EEG_processed_epoched_l1_break_epoch' num2str(k) '.mat']);
    l1_break_epoch = break_epoch_l1.l1_break_epoch;
    instruction_epoch_l2 = load (['C:\Users\Pragati\git\Neuroinformatics\Course_Project_Stroke\Right_paralysis_preprocessed\' 'EEG_processed_epoched_l2_inst_epoch' num2str(k) '.mat']);
    l2_inst_epoch = instruction_epoch_l2.l2_inst_epoch;
    imagery_epoch_l2 = load (['C:\Users\Pragati\git\Neuroinformatics\Course_Project_Stroke\Right_paralysis_preprocessed\' 'EEG_processed_epoched_l2_imagery_epoch' num2str(k) '.mat']);
    l2_imagery_epoch = imagery_epoch_l2.l2_imagery_epoch;
    break_epoch_l2 = load (['C:\Users\Pragati\git\Neuroinformatics\Course_Project_Stroke\Right_paralysis_preprocessed\' 'EEG_processed_epoched_l2_break_epoch' num2str(k) '.mat']);
    l2_break_epoch = break_epoch_l2.l2_break_epoch;
    

    c2plot = 15;

    signal = l1_imagery_epoch;
    baseline = l1_break_epoch;
    
    min_freq =  2;
    max_freq = 80;
    num_frex = 30;

    
    % define wavelet parameters
    time = -1:1/500:1;
    frex = logspace(log10(min_freq),log10(max_freq),num_frex); % 30 equally spacing in the log frequency 
    s = logspace(log10(3),log10(10),num_frex)./(2*pi*frex);
    
%     s = (N/(2*pi*frex))^2; 
    
    n_wavelet            = length(time);
    n_data               = size(signal,2)*size(signal,3);
    n_convolution        = n_wavelet+n_data-1;
    n_conv_pow2          = pow2(nextpow2(n_convolution));
    half_of_wavelet_size = (n_wavelet-1)/2;
    
    
    base = size(baseline,2)*size(baseline,3);
    
    stimuli = squeeze(signal(c2plot,:,:));
    stimuli_baseline = squeeze(baseline(c2plot,:,:));
    % get FFT of data
    eegfft = fft(reshape(stimuli,1,n_data),n_conv_pow2);
    
    eegfft_baseline = fft(reshape(stimuli_baseline,1,base));
    
    
    baseline_power =  mean(abs(reshape(eegfft_baseline,size(baseline,2),size(baseline,3))).^2,2);
    
    baseline_mean = mean(baseline_power(300:700));
    
    % initialize
    eegpower = zeros(num_frex,size(signal ,2)); % frequencies X time X trials
    
    % eegpower_baseline = zeros(num_frex,size(baseline ,2));
    
    % baseidx = dsearchn(EEG.times',[-500 -200]'); %baseline epoch for the task 
    
    % loop through frequencies and compute synchronization
    for fi=1:num_frex
        
        wavelet = wav_amp * fft( sqrt(1/(s(fi)*sqrt(pi))) * exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*(s(fi)^2))) , n_conv_pow2 ); % freqquency of the convolution morlet wavelet 
        % convolution
        eegconv = ifft(wavelet.*eegfft);%element wise freqeuncy domain multuplication 
        eegconv = eegconv(1:n_convolution);
        eegconv = eegconv(half_of_wavelet_size+1:end-half_of_wavelet_size);
        
        
        temppower = mean(abs(reshape(eegconv,size(stimuli,1),size(stimuli,2))).^2,2);
        eegpower(fi,:) = 10*log10(temppower./baseline_mean); 
        %     eegpower(fi,:) = 10*log10(temppower);
    end
    
%     figure
%     subplot(121)
%     contourf(1:size(stimuli,1),frex,eegpower,40,'linecolor','none')
%     set(gca,'clim',[-20 4],'xlim',[0 1000],'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
%     colorbar;
%     title('Logarithmic frequency scaling')
%     
%     subplot(122)
%     contourf(1:size(stimuli,1),frex,eegpower,40,'linecolor','none')
%     colorbar;
%     set(gca,'clim',[-20 4],'xlim',[0 1000])
%     title('Linear frequency scaling')

    if l == 1
        grand_avg_l1 = eegpower ;
    else
        grand_avg_l1 = grand_avg_l1 + eegpower;
    end

    signal = l2_imagery_epoch;
    baseline = l2_break_epoch;
    
    min_freq =  2;
    max_freq = 80;
    num_frex = 30;
    
    % define wavelet parameters
    time = -1:1/500:1;
    frex = logspace(log10(min_freq),log10(max_freq),num_frex); % 30 equally spacing in the log frequency 
    s = logspace(log10(3),log10(10),num_frex)./(2*pi*frex);
    
    n_wavelet            = length(time);
    n_data               = size(signal,2)*size(signal,3);
    n_convolution        = n_wavelet+n_data-1;
    n_conv_pow2          = pow2(nextpow2(n_convolution));
    half_of_wavelet_size = (n_wavelet-1)/2;
    
    
    base = size(baseline,2)*size(baseline,3);
    
    stimuli = squeeze(signal(c2plot,:,:));
    stimuli_baseline = squeeze(baseline(c2plot,:,:));
    % get FFT of data
    eegfft = fft(reshape(stimuli,1,n_data),n_conv_pow2);
    
    eegfft_baseline = fft(reshape(stimuli_baseline,1,base));
    
    
    baseline_power =  mean(abs(reshape(eegfft_baseline,size(baseline,2),size(baseline,3))).^2,2);
    
    baseline_mean = mean(baseline_power(300:700));
    
    % initialize
    eegpower = zeros(num_frex,size(signal ,2)); % frequencies X time X trials
    
    % eegpower_baseline = zeros(num_frex,size(baseline ,2));
    
    % baseidx = dsearchn(EEG.times',[-500 -200]'); %baseline epoch for the task 
    
    % loop through frequencies and compute synchronization
    for fi=1:num_frex
        
        wavelet = wav_amp * fft( sqrt(1/(s(fi)*sqrt(pi))) * exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*(s(fi)^2))) , n_conv_pow2 ); % freqquency of the convolution morlet wavelet 
        % convolution
        eegconv = ifft(wavelet.*eegfft);%element wise freqeuncy domain multuplication 
        eegconv = eegconv(1:n_convolution);
        eegconv = eegconv(half_of_wavelet_size+1:end-half_of_wavelet_size);
        
        
        temppower = mean(abs(reshape(eegconv,size(stimuli,1),size(stimuli,2))).^2,2);
        eegpower(fi,:) = 10*log10(temppower./baseline_mean); 
        %     eegpower(fi,:) = 10*log10(temppower);
    end
    
%     figure
%     subplot(121)
%     contourf(1:size(stimuli,1),frex,eegpower,40,'linecolor','none')
%     set(gca,'clim',[-20 4],'xlim',[0 1000],'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
%     colorbar;
%     title('Logarithmic frequency scaling')
%     
%     subplot(122)
%     contourf(1:size(stimuli,1),frex,eegpower,40,'linecolor','none')
%     colorbar;
%     set(gca,'clim',[-20 4],'xlim',[0 1000])
%     title('Linear frequency scaling')

    if l == 1
        grand_avg_l2 = eegpower / length(participant_id);
    else
        grand_avg_l2 = grand_avg_l2 + eegpower / length(participant_id);
    end
    
end


figure
subplot(121)
contourf(1:size(stimuli,1),frex,grand_avg_l1,40,'linecolor','none')
set(gca,'clim',[-20 4],'xlim',[0 2000],'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colorbar;
title('Average Logarithmic frequency scaling for label 1 across all participants')

subplot(122)
contourf(1:size(stimuli,1),frex,grand_avg_l1,40,'linecolor','none')
colorbar;
set(gca,'clim',[-20 4],'xlim',[0 2000])
title('Average Linear frequency scaling for label 1 across all participants')


figure
subplot(121)
contourf(1:size(stimuli,1),frex,grand_avg_l2,40,'linecolor','none')
set(gca,'clim',[-20 4],'xlim',[0 2000],'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colorbar;
title('Average Logarithmic frequency scaling for label 2 across all participants')

subplot(122)
contourf(1:size(stimuli,1),frex,grand_avg_l2,40,'linecolor','none')
colorbar;
set(gca,'clim',[-20 4],'xlim',[0 2000])
title('Average Linear frequency scaling for label 2 across all participants')


%% STFT with window size of 128 
l1_imag_avg = zeros(1,2000);
l2_imag_avg = zeros(1,2000);

chann = 15;
for l=1:1:length(participant_id_left)
    k = participant_id_left(l);
    
    imagery_epoch_l1= load (['C:\Users\Pragati\git\Neuroinformatics\Course_Project_Stroke\Left_paralysis_preprocessed\' 'EEG_processed_epoched_l1_imagery_epoch' num2str(k) '.mat']);
    l1_imagery_epoch = imagery_epoch_l1.l1_imagery_epoch;
   
    imagery_epoch_l2 = load (['C:\Users\Pragati\git\Neuroinformatics\Course_Project_Stroke\Left_paralysis_preprocessed\' 'EEG_processed_epoched_l2_imagery_epoch' num2str(k) '.mat']);
    l2_imagery_epoch = imagery_epoch_l2.l2_imagery_epoch;
    
   
%     % Plots of spectotopo for each of the 3 stages 
%     figure;
%     subplot(2,1,1);
%     spectopo(l1_inst_epoch,1000,EEG.srate); 
% 
%     subplot(2,1,2);
%     spectopo(l2_inst_epoch,1000,EEG.srate); 
% 
%     figure;
%     subplot(2,1,1);
%     spectopo(l1_imagery_epoch,2000,EEG.srate); 
% 
%     subplot(2,1,2);
%     spectopo(l2_imagery_epoch,2000,EEG.srate); 
% 
%     figure;
%     subplot(2,1,1);
%     spectopo(l1_break_epoch,1000,EEG.srate); 
% 
%     subplot(2,1,2);
%     spectopo(l2_break_epoch,1000,EEG.srate); 
    currsignal =  l1_imagery_epoch;

    

    channel_index = chann;

    l1_imag_avg = l1_imag_avg + squeeze(mean(currsignal(channel_index,:,:),3)) / length(participant_id);

%     figure;
%     plot(1: x_axis_limit,squeeze(mean(currsignal(channel_index,:,:),3)))

    currsignal =  l2_imagery_epoch;
    channel_index = chann;
    x_axis_limit = size(currsignal,2); % in ms

    l2_imag_avg = l2_imag_avg + squeeze(mean(currsignal(channel_index,:,:),3)) / length(participant_id);

%      figure;
%      plot(1: x_axis_limit,squeeze(mean(currsignal(channel_index,:,:),3)))
end


c2plot = chann;% ChannelC3
window = 128;%window_size
fs = 500;%sampling rate

signal = l1_imag_avg;
t = (1:size(signal, 2)) / fs; 
figure;
subplot(2,1,1);
spectrogram(signal, hann(window), window/2, window*2, fs, 'yaxis');
title(['Average Trial Short-Time Fourier Transform (STFT) for  for channel '  num2str(c2plot) ' label-1 Imagery']);
colorbar;
caxis([-15 5])
ylim([0, 30]);
xlabel('Time (s)');
ylabel('Frequency (Hz)');

signal = l2_imag_avg;
t = (1:size(signal, 2)) / fs;

subplot(2,1,2);
spectrogram(signal, hann(window), window/2, window*2, fs, 'yaxis');
title(['Average Trial Short-Time Fourier Transform (STFT) for channel '  num2str(c2plot) ' label-2 Imagery']);
colorbar; 
caxis([-15 5])
ylim([0, 30]); 
xlabel('Time (s)');
ylabel('Frequency (Hz)');


%% ERS/ ERD with 4th order butterworth filter. 


l1_imag_avg = zeros(1,2000);
l2_imag_avg = zeros(1,2000);

chann = 15;
for l=1:1:length(participant_id_right)
    k = participant_id_right(l);
    
    imagery_epoch_l1= load (['C:\Users\Pragati\git\Neuroinformatics\Course_Project_Stroke\Right_paralysis_preprocessed\' 'EEG_processed_epoched_l1_imagery_epoch' num2str(k) '.mat']);
    l1_imagery_epoch = imagery_epoch_l1.l1_imagery_epoch;
   
    imagery_epoch_l2 = load (['C:\Users\Pragati\git\Neuroinformatics\Course_Project_Stroke\Right_paralysis_preprocessed\' 'EEG_processed_epoched_l2_imagery_epoch' num2str(k) '.mat']);
    l2_imagery_epoch = imagery_epoch_l2.l2_imagery_epoch;
    
   
%     % Plots of spectotopo for each of the 3 stages 
%     figure;
%     subplot(2,1,1);
%     spectopo(l1_inst_epoch,1000,EEG.srate); 
% 
%     subplot(2,1,2);
%     spectopo(l2_inst_epoch,1000,EEG.srate); 
% 
%     figure;
%     subplot(2,1,1);
%     spectopo(l1_imagery_epoch,2000,EEG.srate); 
% 
%     subplot(2,1,2);
%     spectopo(l2_imagery_epoch,2000,EEG.srate); 
% 
%     figure;
%     subplot(2,1,1);
%     spectopo(l1_break_epoch,1000,EEG.srate); 
% 
%     subplot(2,1,2);
%     spectopo(l2_break_epoch,1000,EEG.srate); 
    currsignal =  l1_imagery_epoch;

    

    channel_index = chann;

    l1_imag_avg = l1_imag_avg + squeeze(mean(currsignal(channel_index,:,:),3)) / length(participant_id);

%     figure;
%     plot(1: x_axis_limit,squeeze(mean(currsignal(channel_index,:,:),3)))

    currsignal =  l2_imagery_epoch;
    channel_index = chann;
    x_axis_limit = size(currsignal,2); % in ms

    l2_imag_avg = l2_imag_avg + squeeze(mean(currsignal(channel_index,:,:),3)) / length(participant_id);

%      figure;
%      plot(1: x_axis_limit,squeeze(mean(currsignal(channel_index,:,:),3)))
end


ord = 4;
variable_names = { 'l1_imag_avg', 'l2_imag_avg'};   
figure;
for label = 1:1:2
    currsignal = eval(variable_names{label});
    
    channel_index = chan;
    x_axis_limit = size(currsignal,2); %

    % Design the bandpass filter
    fs = 500;
    fpass = [8, 30]; % Passband frequencies in Hz
    order = ord; % Filter order
    [b, a] = butter(order, fpass/(fs/2), 'bandpass'); % Design Butterworth bandpass filter coefficients

    filtered_signal = filter(b, a, currsignal);
%     plot(abs(fft(filtered_signal)));
    smoothness_coefficient = 200;
    smoothed_signal = smoothdata(filtered_signal, 'movmean', smoothness_coefficient);
    
    basemean = mean(smoothed_signal);

    smoothed_signal = smoothed_signal-basemean;
    smoothed_signal=(smoothed_signal/basemean)*100;
 
    plot(1: x_axis_limit,smoothed_signal);
    hold on;

%      figure;
%      for chan= 1:26
%          plot(1: x_axis_limit,squeeze(mean(currsignal(chan,:,:),3)))
%          hold on; 
%       end
    
end


    
    temppower = mean(abs(reshape(eegconv,size(stimuli,1),size(stimuli,2))).^2,2);
    eegpower(fi,:) = 10*log10(temppower./baseline_mean); 
    %     eegpower(fi,:) = 10*log10(temppower);
end

figure
subplot(121)
contourf(1:size(stimuli,1),frex,eegpower,40,'linecolor','none')
set(gca,'clim',[-20 4],'xlim',[0 1000],'yscale','log','ytick',logspace(log10(min_freq),log10(max_freq),6),'yticklabel',round(logspace(log10(min_freq),log10(max_freq),6)*10)/10)
colorbar;
title('Logarithmic frequency scaling')

subplot(122)
contourf(1:size(stimuli,1),frex,eegpower,40,'linecolor','none')
colorbar;
set(gca,'clim',[-20 4],'xlim',[0 1000])
title('Linear frequency scaling')
%% To plot ERS/ERD
% C3-14, C4- 15. 

chan = 15;
variable_names = { 'l1_imagery_epoch', 'l2_imagery_epoch'};   
figure;
for label = 1:1:2
    currsignal = eval(variable_names{label});
    
    channel_index = chan;
    x_axis_limit = size(currsignal,2); %

    % Design the bandpass filter
    fs = 500;
    fpass = [8, 30]; % Passband frequencies in Hz
    order = 1; % Filter order
    [b, a] = butter(order, fpass/(fs/2), 'bandpass'); % Design Butterworth bandpass filter coefficients

    filtered_signal = filter(b, a, mean(currsignal(chan,:,:),3));
    plot(abs(fft(filtered_signal)));
    smoothness_coefficient = 200;
    smoothed_signal = smoothdata(filtered_signal, 'movmean', smoothness_coefficient);
    
    basemean = mean(smoothed_signal);

    smoothed_signal = smoothed_signal-basemean;
    smoothed_signal=(smoothed_signal/basemean)*100;
 
    plot(1: x_axis_limit,smoothed_signal);
    hold on;

%      figure;
%      for chan= 1:26
%          plot(1: x_axis_limit,squeeze(mean(currsignal(chan,:,:),3)))
%          hold on; 
%       end
   
end

%% Performed permutation test to check if there is any significant difference between L1 & L2 and Left and Right Hemiplegia patients. 

k = participant_id_left(l);

imagery_epoch_l1= load (['C:\Users\Pragati\git\Neuroinformatics\Course_Project_Stroke\Left_paralysis_preprocessed\' 'EEG_processed_epoched_l1_imagery_epoch' num2str(k) '.mat']);
l1_imagery_epoch = imagery_epoch_l1.l1_imagery_epoch;


imagery_epoch_l2 = load (['C:\Users\Pragati\git\Neuroinformatics\Course_Project_Stroke\Left_paralysis_preprocessed\' 'EEG_processed_epoched_l2_imagery_epoch' num2str(k) '.mat']);
l2_imagery_epoch = imagery_epoch_l2.l2_imagery_epoch;

data1 = l1_imagery_epoch;
data2 = l2_imagery_epoch;

observed_statistic = mean(data1(:)) - mean(data2(:));

min_epochs = min(size(data1, 3), size(data2, 3));

num_permutations = min_epochs;

permuted_statistics = zeros(num_permutations, 1);

for i = 1:num_permutations
    shuffled_labels = randperm(min_epochs);
    
    shuffled_data1 = data1(:, :, shuffled_labels);
    shuffled_data2 = data2(:, :, shuffled_labels);
    
    permuted_statistics(i) = mean(shuffled_data1(:)) - mean(shuffled_data2(:));
end

p_value = sum(abs(permuted_statistics) >= abs(observed_statistic)) / num_permutations;

disp(['Observed Test Statistic: ', num2str(observed_statistic)]);