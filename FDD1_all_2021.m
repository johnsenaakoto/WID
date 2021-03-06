function hrf=FDD1_all_2021(d,s,dt,WL,TL,LPF)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%water level deconvolution                       %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hrf=FDD_all_2020(d,s,dt,WL,TL,LPF)
% d = radial components or data to deconvolve (i.e. for PP the entire wavetrain)
% s=source funciton (vertical in RF work PP with prevents stuff zeroed out for PP
% dt= sample rate
% WL > 0.101 does curved prewhitening
% WL < 0.101 does a flat line water level deconvolution
% TL=taper length [5 5]
% LPF=frequency for the low pass filter
% meantype='a' for arithmeric mean  'g' for geometric mean (does not work)
%

[rows_d,cols_d]=size(d);
% d=GFILT_data_2021(d,LPF,[20 20],dt,'l',30);
% s=GFILT_data_2021(s,LPF,[20 20],dt,'l',30);

D=fft(d);
S=fft(s);

De=S.*conj(S);
N=D.*conj(S);

DeWL=WL*(ones(rows_d,1)*max(De));
DW=max(De,DeWL);

hrf=real(ifft(N./DW));       %%%%inverse fourier transform%%%%
vrf=real(ifft(De./DW));


hrf = GFILT_data_2021(hrf,LPF,[TL],dt,'l',30);  % turn off and turn off line 47 for one test turn on for other test
% whos
vrf = GFILT_data_2021(vrf,LPF,[TL],dt,'l',30);

nrf=ones(rows_d,1)*(1./max(abs(vrf)));

hrf=nrf.*hrf;

end


