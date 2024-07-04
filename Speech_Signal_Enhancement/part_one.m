
x = load('Sayyam Palrecha - speech_with_beeps.txt');
plotting(x, 1, 'Time', 'Frequency', 'Spectrogram of speech with beeps');
%soundsc(Sayyam_Palrecha___speech_with_beeps);

%STFT algorithm:

% Parameters for computing STFT
fs = 8000; 
winlen = 160; 
overlap = winlen/2;
nfft = winlen; 

% Compute STFT using Hamming window
[S, f, t] = spectrogram(x, hamming(winlen), overlap, nfft, fs);

S_unity = S./abs(S);

% Convert magnitude to dB scale
S_dB = 10*log10(abs(S));

% Plot the spectrogram
specplot(S_dB, t, f, 2, 'Time (s)', 'Frequency (Hz)', 'Spectrogram of STFT');

%peak detection algorithm:

avg_value = mean(S_dB);
sigma = std(S_dB);

%tuning parameter alpha

window_size = winlen;
alpha = ones(1, size(S_dB, 2)); 
smooth_S_dB = adaptive_smoothing(S_dB, window_size);

% Compute alpha based on the difference between original and smoothed signal
for i = 1:size(S_dB, 2)
    alpha(1, i) = abs(mean(S_dB(:, i) - smooth_S_dB(:, i)))./max(abs(mean(S_dB(:, i) - smooth_S_dB(:, i))))/2;
end

threshold = avg_value + alpha.*sigma; % minimum requirement for flagging a sample as a part of a beep

data_set = []; % extract the beep data which includes time and frequency indices along with the energy levels
[row, col, val, data_set] = detect_beep(S_dB, data_set, threshold);

maxspace = floor(sigma./alpha); % maximum time samples spacing allowed between two consecutive detected samples
minlen = (maxspace).^2; % minimum length of a beep signal to be detected
% maximum limit for the count of samples allowed in one round in the interval of an expected beep
maxlim = floor(max(histcounts(row, 'BinMethod','integer'))./sigma./alpha);

%retrieving the speech signal which was detected as beep
for i = 1:length(f)
    if ~isempty(find(row == i, 1))
        t_srt = find(row == i, 1);
        t_end = find(row == i, 1) + sum(row == i) - 1;
        count = 0;
        for j = t_srt:t_end - 1
            if col(j+1) - col(j) < maxspace(col(j))
                count = count + 1;
            elseif count < minlen(col(j))
                for k = 1:count+1
                        val(j - k + 1) = 0;
                end
                count = 0;
            end
            if count >= maxlim(col(j)) && val(j) < (median(val))./alpha(col(j))/2
                count = 0;
            end
        end
    end
end

%updating the data_set
data_set(:, 3) = val;

beep_data = zeros(size(S_dB));

% assigning values from 'val' array at corresponding positions in the 'beep_data' array
beep_data(sub2ind(size(beep_data), row, col)) = val; 

specplot(beep_data, t, f, 3, 'Time (s)', 'Frequency (Hz)', 'Spectrogram of beep signal')

% adjusted gain based on the statistical information of the beep_data, so that it can
% be adjusted for an arbitary signal containing beeps.
Gain = sigma./alpha;
speech_signal = S_dB - Gain.*beep_data;

%algorithm for reconstruction of the speech signal here:-

speech_abs = 10.^(speech_signal/10);
speech_comp = speech_abs.*S_unity;

y = istft(speech_comp, winlen, overlap, nfft);

%plotting spectrofram of speech signal
plotting(y, 4, 'Time (s)', 'Frequency (Hz)', 'Spectrogram of speech signal');

soundsc(y);
filename = 'speech_without_beeps.wav';  
audiowrite(filename, y, fs, 'BitsPerSample', 24);

% ISTFT function
function x = istft(S, window, overlap, nfft)
    xlen = (size(S,2)-1)*overlap + nfft;
    x = zeros(xlen, 1);
    for n = 1:size(S,2)
        start_idx = (n-1)*overlap + 1;
        end_idx = start_idx + nfft - 1;
        x_seg = real(ifft(S(:,n), nfft));
        x(start_idx:end_idx) = x(start_idx:end_idx) + x_seg .* window;
    end
    x = (x/sum(window))'; %scaling
end

function plotting(y, index, xaxis, yaxis, Title)
    subplot(2,2,index);
    specgram(y);
    xlabel(xaxis);
    ylabel(yaxis);
    title(Title);
end

function specplot(energy, time, frequency, index, xaxis, yaxis, Title)
    subplot(2,2,index);
    imagesc(time, frequency, energy); 
    axis xy;
    colormap(jet); 
    xlabel(xaxis);
    ylabel(yaxis);
    title(Title);
end

%algorithm for detecting beeps
function [row, col, val, data_set] = detect_beep(S_dB, data_set, threshold)
    for i = 1:size(S_dB, 1)
        for j = 1:size(S_dB, 2)
            if S_dB(i,j) > threshold 
                data_set = [data_set; [i, j, S_dB(i, j)]];
            end
        end
    end

    row = data_set(:, 1);
    col = data_set(:, 2);
    val = data_set(:, 3);
end

function smoothed_signal = adaptive_smoothing(signal, window_size) 
    smoothed_signal = zeros(size(signal));
    for i = 1:length(signal)
        start_idx = max(1, i - floor(window_size / 2));
        end_idx = min(length(signal), i + floor(window_size / 2));
        smoothed_signal(i) = mean(signal(start_idx:end_idx));
    end
end
