def WID_FID(sig, src, num_it, dt,LPF,dzz,dtt,len_data, WL)
    #WID_FID takes in signal and source fuction pairs from a
    #CMP gather and performs wavefield iterative deconvolution on them to
    #produce a stacked SS receiver function

    #The input arguments are:
    #signal=array of signal functions
    #source=array of source functions
    #num_it=number of iterations for deconvolution
    #dt=sampling rate in seconds/sample
    #LPF=low pass frequency for filter
    #dzz=depths for signal functions obtained from a velocity model
    #dtt=travel times for signal functions obtained from a velocity model
    #phase_parameters= phase info from deconvolution
    #WL =Water level for deconvolution

    [nrow, ncol] = 
