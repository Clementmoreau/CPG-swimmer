function Adat_filtered = filter_high_mode(Adat,Nfilt)

% Adat : data array made of columns of length n
% Nfilt : number of modes kept

% Fourier transform
Adatf = fft(Adat);

% Prepare the filtered data
Adatf_filt = zeros(size(Adat));

% filter high modes
Adatf_filt(1:Nfilt+1,:) = Adatf(1:Nfilt+1,:);
Adatf_filt(end-Nfilt+1:end,:) = Adatf(end-Nfilt+1:end,:);

% inverse fourier transform
Adat_filtered = ifft(Adatf_filt,'symmetric');

end
