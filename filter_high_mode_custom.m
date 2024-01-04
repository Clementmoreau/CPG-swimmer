function Adat_filtered = filter_high_mode_custom(Adat,Nfilt,N)

Adatf = fft(Adat);
Adatf_filt = zeros(size(Adat));

% filter high modes
Adatf_filt(1:Nfilt+1,:) = Adatf(1:Nfilt+1,:);
Adatf_filt(end-Nfilt+1:end,:) = Adatf(end-Nfilt+1:end,:);

% keep also the 2 highest modes
Adatf_filt(N/2:N/2+2,:) = Adatf(N/2:N/2+2,:);

% inverse fourier transform
Adat_filtered = ifft(Adatf_filt,'symmetric');

end
