clear all;
path0=pwd;
path2=[ '../../Nonequilibrium'];
addpath(path2);

for FLAG=1:4
    FLAG
    indexN=[1:31 50:80];
    
    N=62;
    Tau=2;
    NSUB=176;
    
    % Parameters of the data
    TR=1;  % Repetition Time (seconds)
    % Bandpass filter settings
    fnq=1/(2*TR);                 % Nyquist frequency
    flp = 0.008;                    % lowpass frequency of filter (Hz)
    fhi = 0.08;                    % highpass
    Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
    k=2;                          % 2nd order butterworth filter
    [bfilt,afilt]=butter(k,Wn);   % construct the filter
    
    
    if FLAG==1
        load hcp7t_rfMRI_REST2_AP_dbs62.mat;
    elseif FLAG==2
        load hcp7t_rfMRI_REST1_PA_dbs62.mat;
    elseif FLAG==3
        load hcp7t_rfMRI_REST3_PA_dbs62.mat;
    elseif FLAG==4
        load hcp7t_rfMRI_REST4_AP_dbs62.mat;
    end
    
    FowRev_R=zeros(1,NSUB);
    FowRev_M=zeros(1,NSUB);
    
    for sub=1:NSUB  % over subjects
        ts2=subject{sub}.dbs80ts;
        ts=ts2(indexN,:);
        clear signal_filt;
        for seed=1:N
            ts(seed,:)=detrend(ts(seed,:)-nanmean(ts(seed,:)));    %filtering
            signal_filt(seed,:) =filtfilt(bfilt,afilt,ts(seed,:));
        end
        ts=signal_filt(:,10:end-10);
        Tm=size(ts,2);
        FCtf=corr(ts(:,1:Tm-Tau)',ts(:,1+Tau:Tm)');       %% Core...FC tau foward
        FCtr=corr(ts(:,Tm:-1:Tau+1)',ts(:,Tm-Tau:-1:1)'); %% FC tau reversal
        Itauf=-0.5*log(1-FCtf.*FCtf);  %% Mutual information...
        Itaur=-0.5*log(1-FCtr.*FCtr);
        Reference=((Itauf(:)-Itaur(:)).^2)';
        index=find(Reference>quantile(Reference,0.0));
        FowRev_R(sub)=nanmean(Reference(index));
        Tenet2_R(sub,:,:)=abs(Itauf-Itaur);
    end
    
    if FLAG==1
        load hcp7t_tfMRI_MOVIE1_AP_dbs62.mat;
    elseif FLAG==2
        load hcp7t_tfMRI_MOVIE2_PA_dbs62.mat;
    elseif FLAG==3
        load hcp7t_tfMRI_MOVIE3_PA_dbs62.mat;
    elseif FLAG==4
        load hcp7t_tfMRI_MOVIE4_AP_dbs62.mat;
    end
    
    for sub=1:NSUB
        ts2=subject{sub}.dbs80ts;
        ts=ts2(indexN,:);
        clear signal_filt;
        for seed=1:N
            ts(seed,:)=detrend(ts(seed,:)-nanmean(ts(seed,:)));
            signal_filt(seed,:) =filtfilt(bfilt,afilt,ts(seed,:));
        end
        ts=signal_filt(:,10:end-10);
        Tm=size(ts,2);
        FCtf=corr(ts(:,1:Tm-Tau)',ts(:,1+Tau:Tm)');
        FCtr=corr(ts(:,Tm:-1:Tau+1)',ts(:,Tm-Tau:-1:1)');
        Itauf=-0.5*log(1-FCtf.*FCtf);
        Itaur=-0.5*log(1-FCtr.*FCtr);
        Reference=((Itauf(:)-Itaur(:)).^2)';
        index=find(Reference>quantile(Reference,0.0));
        FowRev_M(sub)=nanmean(Reference(index));
        Tenet2_M(sub,:,:)=abs(Itauf-Itaur);
    end
    
    THRLOW=0;
    THRHIGH=100;
    FowRev_R=rmoutliers(FowRev_R,'percentiles',[THRLOW THRHIGH]);
    FowRev_M=rmoutliers(FowRev_M,'percentiles',[THRLOW THRHIGH]);
    
    figure(FLAG);
    violinplot([FowRev_R' FowRev_M']);
    ranksum(FowRev_R,FowRev_M)
    
    if FLAG==1
        save results_MOVIE1_REST2_62.mat FowRev_R FowRev_M;
    elseif FLAG==2
        save results_MOVIE2_REST1_62.mat FowRev_R FowRev_M;
    elseif FLAG==3
        save results_MOVIE3_REST3_62.mat FowRev_R FowRev_M;
    elseif FLAG==4
        save results_MOVIE4_REST4_62.mat FowRev_R FowRev_M;
    end
    
end
