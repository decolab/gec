#!/bin/bash
#SBATCH --qos=vip
#SBATCH --job-name=Tenet3
#SBATCH --mail-type=END
#SBATCH --mail-user=
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=2
#SBATCH --array=1-180
#SBATCH --output=Tenet3%A_%a.out
#SBATCH --error=Tenet3%A_%a.err

#Load Matlab 2017a module
ml MATLAB

matlab -nojvm -nodisplay<<-EOF

s=str2num(getenv('SLURM_ARRAY_TASK_ID'))

load /home/gdeco/NonEquilibrium/hcp7t_tfMRI_MOVIE1_AP_dbs62.mat;

load /home/gdeco/NonEquilibrium/SC_dbs80HARDIFULL.mat;

indexregion=[1:31 50:80];

C=SC_dbs80HARDI;
C=C/max(max(C))*0.2;
C=C(indexregion,indexregion);

FLAGREV=1;

NSUB=181;
NSIM=100;
N=62;
Tau=2;
ThRef=0.75;

% Parameters of the data
TR=1;  % Repetition Time (seconds)

% Bandpass filter settings
fnq=1/(2*TR);                 % Nyquist frequency
flp = 0.008;                    % lowpass frequency of filter (Hz)
fhi = 0.08;                    % highpass
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);   % construct the filter

Isubdiag = find(tril(ones(N),-1));

%%%%%%%%%%%%%%

nsub=1;
for sub=1:NSUB
    if isstruct(subject{sub})
        if size(subject{sub}.dbs80ts,2)==size(subject{1}.dbs80ts,2)
            subject{nsub}=subject{sub};
            nsub=nsub+1;
        end
    end
end

NSUB=nsub-1;

%%%%%%%%%%%%%%

ts=subject{s}.dbs80ts;
ts=ts(indexregion,:);
Tmax=size(ts,2);
signal_filt2=zeros(N,Tmax);

for seed=1:N
    ts(seed,:)=detrend(ts(seed,:)-nanmean(ts(seed,:)));
    signal_filt2(seed,:)=filtfilt(bfilt,afilt,ts(seed,:));
end
signal_filt=signal_filt2(:,20:end-20);
Tm=size(signal_filt,2);

FCemp=corrcoef(signal_filt');
Ifcemp=-0.5*log(1-FCemp.*FCemp);


%% FC(tau)
FCtaufemp=corr(signal_filt(:,1:Tm-Tau)',signal_filt(:,1+Tau:Tm)');
FCtauremp=corr(signal_filt(:,Tm:-1:Tau+1)',signal_filt(:,Tm-Tau:-1:1)');

Itaufemp=-0.5*log(1-FCtaufemp.*FCtaufemp);
Itauremp=-0.5*log(1-FCtauremp.*FCtauremp);

f_diff=0.018*ones(1,N);

% Parameters HOPF
a=-0.02*ones(N,2);
omega = repmat(2*pi*f_diff',1,2); omega(:,1) = -omega(:,1);
dt=0.1*TR/2;
sig=0.01;
dsig = sqrt(dt)*sig;

%%%%%%%%%%%%
%% Optimize
%%

FCsim2=zeros(NSIM,N,N);
FCtaufsim2=zeros(NSIM,N,N);
FCtaursim2=zeros(NSIM,N,N);
Ifcs1=zeros(NSIM,N,N);
Itaufs1=zeros(NSIM,N,N);
Itaurs1=zeros(NSIM,N,N);

we=1;
Cnew=zeros(N,N);
for iter=1:3000
    iter
    wC = we*Cnew;
    sumC = repmat(sum(wC,2),1,2); % for sum Cij*xj
    for sub=1:NSIM
        xs=zeros(Tmax,N);
        z = 0.1*ones(N,2); % --> x = z(:,1), y = z(:,2)
        nn=0;
        % discard first 3000 time steps
        for t=0:dt:3000
            suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
            zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
            z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(N,2);
        end
        % actual modeling (x=BOLD signal (Interpretation), y some other oscillation)
        for t=0:dt:((Tmax-1)*TR)
            suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
            zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
            z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(N,2);
            if abs(mod(t,TR))<0.01
                nn=nn+1;
                xs(nn,:)=z(:,1)';
            end
        end
        
        %%%%
        BOLD=xs';
        signal_filt22=zeros(N,nn);
        for seed=1:N
            BOLD(seed,:)=demean(detrend(BOLD(seed,:)));
            signal_filt22(seed,:) =filtfilt(bfilt,afilt,BOLD(seed,:));
        end
        signal_filt=signal_filt22(:,20:end-20);
        Tm=size(signal_filt,2);
        FCsim2(sub,:,:)=corrcoef(signal_filt');
        Ifcs1(sub,:,:)=-0.5*log(1-squeeze(FCsim2(sub,:,:)).*squeeze(FCsim2(sub,:,:)));
        
        %% FC(tau)
        
        FCtaufsim2(sub,:,:)=corr(signal_filt(:,1:Tm-Tau)',signal_filt(:,1+Tau:Tm)');
        FCtaursim2(sub,:,:)=corr(signal_filt(:,Tm:-1:Tau+1)',signal_filt(:,Tm-Tau:-1:1)');
        
        Itaufs1(sub,:,:)=-0.5*log(1-squeeze(FCtaufsim2(sub,:,:)).*squeeze(FCtaufsim2(sub,:,:)));
        Itaurs1(sub,:,:)=-0.5*log(1-squeeze(FCtaursim2(sub,:,:)).*squeeze(FCtaursim2(sub,:,:)));
        
        %% Reversibility TENET2
        
        Itauf=-0.5*log(1-squeeze(FCtaufsim2(sub,:,:)).*squeeze(FCtaufsim2(sub,:,:)));
        Itaur=-0.5*log(1-squeeze(FCtaursim2(sub,:,:)).*squeeze(FCtaursim2(sub,:,:)));
        Reference=((Itauf(:)-Itaur(:)).^2)';
        index=find(Reference>quantile(Reference,ThRef));
        FowRevsim2(sub)=sqrt(mean(Reference(index)));
    end
    
    Ifcsim=squeeze(mean(Ifcs1));
    Itaufsim=squeeze(mean(Itaufs1));
    Itaursim=squeeze(mean(Itaurs1));
    
    for i=1:N
        for j=1:N
            if (C(i,j)>0 || j==N-i+1)
                Cnew(i,j)=Cnew(i,j)+0.0005*(Ifcemp(i,j)-Ifcsim(i,j)) ...
                    -FLAGREV*0.0001*(Itaufemp(i,j)-Itaufsim(i,j)) ...
                    +FLAGREV*0.0001*(Itauremp(i,j)-Itaursim(i,j));
                if Cnew(i,j)<0
                    Cnew(i,j)=0;
                end
            end
        end
    end
end

EffectiveConnectivity=Cnew;

save (sprintf('results_eff_hopf_fcrev_movie_MOVIE1sub_%03d.mat',s),'FCemp','EffectiveConnectivity');

EOF
