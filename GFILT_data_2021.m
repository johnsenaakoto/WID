function d=GFILT_data_2021(data,f,TW,dt,flag,pad)
%This program will detrend, taper, then apply a gaussian filter to your
%seismic data.

%PAD FRONT AND BACK WITH ZEROS
%data=[zeros(round(30/dt),1); data; zeros(round(30/dt),1)];

[rows_d,cols_d]=size(data);
if pad>1
    data=[zeros(round(pad/dt),cols_d); data; zeros(round(pad/dt),cols_d)];
end
[rows_d,cols_d]=size(data);


%runs loop for all components
data=detrend(data); %detrends seismogram data
taper=ones(1,rows_d);
ITW=round(TW./dt);

if TW(1)>0.5
    taper(1:ITW(1)+1)=[0:1/ITW(1):1];
    taper(end-ITW(2):end)=[1:-1/ITW(2):0];
    taper2=taper'*ones(1,cols_d);
    data=taper2.*data;
end

a=2*f; %a=2*f is .1 power in a gaussian for gaussian filter formula
Fn=1/2/dt; %nyquist frequency
freq=([0:Fn/(ceil(rows_d/2)-1):Fn]);
%setting up n/2 frequencies between 0 and the nyquist frequencies (this is n frequencies when negative frequencies are included
freq=([freq -fliplr(freq(2:end))]); %complicated way of going from 0 to nyquist frequency, then from -nyquist frequency back to 0
w=2*pi*freq; %angular frequency

if lower(flag(1))=='h'; % high pass filter
    Gw=(1-exp(-1*(w.^2)./(4*a^2)))'; %high pass filter
else
    Gw=(exp(-1*(w.^2)./(4*a^2)))'; %Low pass filter
end

Gw=Gw*ones(1,cols_d);

DATA=fft2(data);
DATA=Gw.*DATA;
d=real(ifft2(DATA));
d = d(1 + round(30/dt):size(d,1) - round(30/dt),:);

end


