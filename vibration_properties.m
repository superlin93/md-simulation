function [attempt_freq, vibration_amp, std_attempt_freq] = vibration_properties(sim_data, show_pics)
% Get the attempt frequency and vibration amplitude
    % Define some stuff (from Matlab example of fft):
    ts = sim_data.time_step; % Sampling period
    fs = 1/ts; %Sampling frequency
    length = sim_data.nr_steps;
    if mod(length, 2) ~= 0 %necessary for one_sided
        length = length - 1;
    end
    one_sided = zeros(sim_data.nr_diffusing, (length/2)+1);
    
%% For each diffusing atom get the speed at each time (derivative of displacement)
%   Then apply the fourier transform
    speed = zeros(sim_data.nr_diffusing, sim_data.nr_steps);
    freq_mean = zeros(sim_data.nr_diffusing, 1);
    atom = 0;
    vib_count = 0;
    for i = sim_data.start_diff_elem:sim_data.end_diff_elem %Loop over diffusing atoms
        atom = atom + 1;
        for time = 2:sim_data.nr_steps
            % speed is in Angstrom/step
            speed(atom, time) = sim_data.displacement(i, time) - sim_data.displacement(i, time-1);
            % Integrate to get vibration amplitude as long as the sign stays the same:
            if sign(speed(atom,time)) == sign(speed(atom,time-1))
                amplitude(vib_count) = amplitude(vib_count) + speed(atom,time);
            else % When the sign changes start the next vibration amplitude
                vib_count = vib_count + 1;
                amplitude(vib_count) = speed(atom,time);
            end
        end
        freq_mean(atom) = meanfreq(speed(atom,:), fs);

        % Fourier transform speed, and get the frequency   
        trans = fft(speed(atom,:));
        two_sided = abs(trans/length);
        one_sided(atom,:) = two_sided(1:length/2 +1);
        one_sided(atom, 2:end-1) = 2*one_sided(atom, 2:end-1); %to get the right amplitude...
    end
    %% Extract the vibration amplitude   
   [mean_vib, vibration_amp] = normfit(amplitude);
   
    %% Extract attempt frequency:
    attempt_freq = mean(freq_mean);
    std_attempt_freq = std(freq_mean);
   
%% Pictures: 
    if show_pics
        figure     % Plot the histogram of vibrational amplitudes
        h = histfit(amplitude, 100); 
        hold on
        %plot([-vibration_amp, vibration_amp], [8000 8000], 'g', 'LineStyle', ':', 'LineWidth', 3)
        %xlim([-2 2])
        h(2).LineWidth = 3;
        title('Histogram of vibrational amplitudes with fitted Gaussian')
        xlabel('Amplitude (Angstrom)')
        ylabel('Occurrence (a.u.)')     
        set(gca, 'YTickLabel','')
        hold off
        
    % Plot the obtained frequency:
        f = fs*(0:(length/2))/length;
        figure
        hold on
        sum_freqs = sum(one_sided);
        %plot(f,sum_freqs)
        % Smooth to  make it look nicer
        smoothed = smooth(sum_freqs, 51);
        plot(f, smoothed, 'LineWidth', 3)
        %title('Frequency spectrum of diffusing element')
        xlabel('Frequency (Hz)')
        ylabel('Occurrence (a.u.)')
        set(gca, 'YTickLabel','')
        % Plot the attempt frequency:
        plot([attempt_freq attempt_freq],0:1, '-r', 'LineWidth', 3)
        % Plot the standard deviations:
        plot([attempt_freq-std_attempt_freq attempt_freq-std_attempt_freq],0:1, ':r', 'LineWidth', 3)       
        plot([attempt_freq+std_attempt_freq attempt_freq+std_attempt_freq],0:1, ':r', 'LineWidth', 3)
        % Limits:
        ylim([0 max(sum_freqs)]) %*1.1])
        xlim([0 2.5E13]) %5eE13 as the maximum should be enough usually
        hold off
    end
end