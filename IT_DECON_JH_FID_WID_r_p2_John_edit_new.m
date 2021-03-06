function [synth_now_stack_z,PP_SS_now_stack_z,phase_info]  = IT_DECON_JH_FID_WID_r_p2_John_edit_new(signal, source, num_it, dt, LPF, dzz, dtt, len_data, phase_info,WL)
%IT_DECON_JH_CC_WID_r_p2 takes in signal and source fuction pairs from a
%CMP gather and performs wavefield iterative deconvolution on them to
%produce a stacked SS receiver function

%The input arguments are:
%signal=array of signal functions
%source=array of source functions
%num_it=number of iterations for deconvolution
%dt=sampling rate in seconds/sample
%LPF=low pass frequency for filter
%dzz=depths for signal functions
%dtt=travel times for signal functions
%phase_parameters= phase info from deconvolution
%WL =Water level for deconvolution
[nrow,ncol]=size(signal);
T=(1/LPF)/dt;
signal_now=signal;
synth_now=zeros(nrow,ncol);
synth_now_stack_z=zeros(1001,1);
phase_info.num_iter= num_it;


source_ac=FDD1_all_2021(source,source,dt,WL,[0 0],LPF);
for n1=1:ncol
    dtt_resampled(:,n1)=interp1(dzz(:,n1),dtt(:,n1),[0:1000]);
    [amp_ac(n1), i_ac(n1)]=max(abs(source_ac(:,n1)));
    sac=source_ac(:,n1);
    sac=flipud(sac(300:(length(sac))));
    source_ac_z(:,n1)=interp1([0:(length(sac)-1)]*dt, sac, dtt_resampled(:,n1));
end

% for n1=1:ncol
%     dtt_resampled(:,n1)=interp1(dzz(:,n1),dtt(:,n1),[0:1000]);
%     source_ac(:,n1)=FDD1_all_2020NSH(source(:,n1),source(:,n1),dt,WL,[0 0],LPF);
%     [amp_ac(n1), i_ac(n1)]=max(abs(source_ac(:,n1)));
%     sac=source_ac(:,n1);
%     sac=flipud(sac(300:(length(sac))));
%     source_ac_z(:,n1)=interp1([0:(length(sac)-1)]*dt, sac, dtt_resampled(:,n1));
% end

source_ac_stack_z=mean(source_ac_z,2);
[amp_ac_stack_z, ~] = max(abs(source_ac_stack_z));

% figure(68)
% plot(source_ac_stack_z,'linewidth',2)
% title('Auto-correlation Depth Stack')
% xlabel('Depth (km)')
% ylabel('Amplitude')
% grid on



for i=1:num_it
    
    sig_cc = FDD1_all_2021(signal_now,source,dt,WL,[0 0],LPF);
    for n2=1:ncol
        DCC=sig_cc(:,n2);
        DCC=flipud(DCC(1:(length(sig_cc(:,n2))-300)));
        sig_cc_z(:,n2)=interp1([0:(length(DCC)-1)]*dt, DCC, dtt_resampled(:,n2));
    end
    
%     for n2=1:ncol
%         sig_cc(:,n2)=FDD1_all_2020NSH(signal_now(:,n2),source(:,n2),dt,WL,[0 0],LPF);
%         DCC=sig_cc(:,n2);
%         DCC=flipud(DCC(1:(length(sig_cc(:,n2))-300)));
%         sig_cc_z(:,n2)=interp1([0:(length(DCC)-1)]*dt, DCC, dtt_resampled(:,n2));
%     end
    
    
    sig_cc_stack_z=mean(sig_cc_z,2);
    [~, i_cc_stack_z] = max(abs(sig_cc_stack_z));
    i_cc_stack_z=i_cc_stack_z-1;
%     figure(69)
%     plot(sig_cc_stack_z,'linewidth',2)
%     title('Cross-correlation Depth Stack')
%     xlabel('Depth (km)')
%     ylabel('Amplitude')
%     grid on

    if i_cc_stack_z<=0
        i_cc_stack_z=1;
    end
    synth_now_stack_z(i_cc_stack_z)=synth_now_stack_z(i_cc_stack_z) + sig_cc_stack_z(i_cc_stack_z)/amp_ac_stack_z;
%     figure(70)
%     plot(synth_now_stack_z, 'linewidth',3)
%     title('Common Stack SdS Depth Function')
%     xlabel('Depth (km)')
%     ylabel('Amplitude')
%     grid on
%     pause

    phase_info.PPFcommon_depth(i)= i_cc_stack_z;
    phase_info.PPFcommon_amp_depth(i)= synth_now_stack_z(i_cc_stack_z);
    
    %scale = 2;
    %time = [1:length(synth_now(:,1))]*dt;
    for n3=1:ncol
        tn=interp1([0:1000],dtt_resampled(:,n3),i_cc_stack_z);
        itn=len_data-round(tn/dt);
        t1=floor(itn-T);
        t2=ceil(itn+T);
        
        if t1<0
            t1=1;
        end
        if t2>length(sig_cc(:,n3))
            t2=length(sig_cc(:,n3));
        end
        
        sig_cc_temp=sig_cc(:,n3);
        
        if sig_cc_stack_z(i_cc_stack_z)>0
            [amp_cc,i_cc] = max((sig_cc_temp(t1:t2)));
        else
            [amp_cc,i_cc] = min((sig_cc_temp(t1:t2)));
        end
        i_cc=i_cc+t1-1;
        i_now=i_cc;
        
        synth_now(i_now,n3)=synth_now(i_now,n3) + amp_cc/2*amp_ac(n3);
%         figure(71)
%                 %plot(time-90, scale*synth_now(:, n3)+n3, 'b', 'linewidth',1)
%                 plot(scale*synth_now(:, n3)+n3, 'b', 'linewidth',1)
%                 title('Individual SdS Time Functions')
%                 xlabel('Time (s)')
%                 ylabel('Number of CMPs')
%                 grid on
%                 hold on
%                 ylim([0 ncol+1])
%                 %axis([100 450 0 (ncol+1)])
        
        sig_synth=convFD(synth_now(:,n3), source(:,n3));
%         figure(72)
%                 %plot(time, scale*signal(:,n3) + n3, 'g','linewidth',3)
%                 plot(scale*signal(:,n3)+n3, 'g','linewidth',3)
%                 hold on
%                 figure(72)
%                 plot(scale*sig_synth + n3, 'k','linewidth',1)
%                 figure(72)
%                 %plot(time, scale*sig_synth + n3,'k','linewidth',1)
%                 %plot(source(:,n3),'k','linewidth',1)
%                 grid on
%                 title('Synthetic Signal and Original Signal')
%                 legend('original signal', 'synthetic signal', 'Location', 'northwest')
%                 xlabel('Time (s)')
%                 ylabel('Number of CMPs')
%                 ylim([0 ncol+1])
%                 %axis([100 450 0 (ncol+1)])
        
        signal_now(:,n3)=signal(:,n3)-sig_synth;
%         figure(73)
%                 %plot(time, scale*signal(:,n3) + n3,'r','linewidth',3)
%                 plot(scale*signal(:,n3)+n3,'r','linewidth',3)
%                 grid on
%                 hold on
%                 %plot(time, scale*signal_now(:,n3) + n3,'k','linewidth',1)
%                 plot(scale*signal_now(:,n3)+n3,'k','linewidth',1)
%                 title('Signal Remaining and Original Signal')
%                 legend('original signal', 'signal remaining', 'Location', 'northwest')
%                 xlabel('Time (s)')
%                 ylabel('Number of CMPs')
%                 %axis([120 450 0 (ncol+1)])
%                 ylim([0 ncol+1])
        
        
        phase_info.S(n3).PPFindividual_time(i)= i_now*dt;   %all relevant phase info is saved into phase_info
        phase_info.S(n3).PPFcommon_time(i)= (itn*dt);
        phase_info.S(n3).PPFindividual_amp(i)= synth_now(i_now,n3);
        %scale = 5;
    end 
%     pause
%     figure(71); clf
%     figure(72); clf
%     figure(73); clf    
end

%close all
for n4=1:ncol
    PP_SS_cc=synth_now(:,n4);
    PP_SS_cc=flipud(PP_SS_cc(1:(length(PP_SS_cc)-300)));
    PP_SS_cc=GFILT_data_2021(PP_SS_cc,1,[0 0],dt,'l', 30);
    PP_SS_cc_z(:,n4)=interp1([0:(length(PP_SS_cc)-1)]*dt, PP_SS_cc, dtt_resampled(:,n4)); %individual SS functions depth converted

%     figure(75)
%     plot((scale*2)*PP_SS_cc_z(:, n4) + n4, 'b', 'linewidth',1)
%     title('Individual SdS Depth Function')
%     xlabel('Depth')
%     ylabel('Number of CMPs')
%     axis([0 1000 0 (ncol+1)])
%     grid on
%     hold on
end


PP_SS_now_stack_z=mean(PP_SS_cc_z,2);
% figure(74)
%     plot(PP_SS_now_stack_z, 'linewidth',3)
%     title('Individual Stack SdS Receiver Function')
%     xlabel('Depth (km)')
%     ylabel('Amplitude')
%     xlim([0 1000])
%     grid on
%     pause
% close all

end


function f0=convFD(f1,f2)
F1=fft(f1);
F2=fft(f2);
f0=real(ifft(F1.*F2));
end