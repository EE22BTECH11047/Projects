
x_noisy = load('Sayyam Palrecha - speech_noisy.txt');

fs = 8000; 
windowLength = 160; 
overlap = windowLength/2; 
nfft = windowLength;

% Compute STFT
[S, f, t] = spectrogram(x_noisy, hamming(windowLength), overlap, nfft, fs);

% PSD estimate
N_frames = size(S,2);
N_freq_bins = size(S,1);
Rxx = zeros(N_frames, N_freq_bins);
alpha = 0.85; %smoothening factor
Rxx(1,:) = (1-alpha)*abs(S(:,1)).^2;
for n = 2:N_frames
    Rxx(n,:) = alpha*Rxx(n-1,:).' + (1-alpha)*abs(S(:,n)).^2;
end

S_squared = abs(S).^2;

%minimum noise power estimation 
r_fraction = 0.001; 
r = round(r_fraction * N_frames);
beta = 1.5; %over-estimation factor 
N_min_squared = zeros(N_frames, N_freq_bins);
for n = 1:N_frames
    for omega_i = 1:N_freq_bins
        start_frame = max(1,n-r);
        Rxx_min = min(Rxx(start_frame:n, omega_i));
        X_squared = S_squared(omega_i, n);
        N_min_squared(n, omega_i) = min(Rxx_min^2, X_squared);
    end
end
N_est_squared=beta*N_min_squared;
%Applying weights
filter = (S_squared - N_est_squared.') ./ max((S_squared+N_est_squared.'), 1e-10);
filter(S_squared < N_est_squared.') = 0;
filtered_S = S .* (filter);

filtered_x = istft(filtered_S, hamming(windowLength), overlap, nfft, fs);

% Plot original and filtered signals
figure;

% Original signal
subplot(2,1,1);
plot(x_noisy);
title('Original Speech Signal');
xlabel('Sample');
ylabel('Amplitude');

% Filtered signal
subplot(2,1,2);
plot(filtered_x);
title('Filtered Speech Signal');
xlabel('Sample');
ylabel('Amplitude');
soundsc(filtered_x, fs);

% ISTFT function
function x = istft(S, window, overlap, nfft, fs)
    xlen = (size(S,2)-1)*overlap + nfft;
    x = zeros(xlen, 1);
    w = window(:); 
    
    for n = 1:size(S,2)
        start_idx = (n-1)*overlap + 1;
        end_idx = start_idx + nfft - 1;
        x_seg = real(ifft(S(:,n), nfft));
        x(start_idx:end_idx) = x(start_idx:end_idx) + x_seg .* w;
    end
    x_scaling = sum(w.^2);
    x = x(1:xlen-nfft+1);
    x = x / x_scaling;
end
