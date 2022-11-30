clear all;

for s=1:180
    load(sprintf('results_eff_hopf_fcrev_movie_MOVIE1sub_%03d.mat',s));
    DataMOVIE1(s,:)=EffectiveConnectivity(:);
end

for s=1:177
    load(sprintf('results_eff_hopf_fcrev_movie_MOVIE4sub_%03d.mat',s));
    DataREST2(s,:)=EffectiveConnectivity(:);
end

task_nsubn={180,177};

trainno=167;
valno=10;

kfold=100;

cl=1:2;
pc=zeros(2,2);

for nfold=1:kfold
    xx=1;
    nfold
    shuffling=randperm(task_nsubn{xx});
    Data=DataMOVIE1;
    Data=Data(shuffling(1:task_nsubn{xx}),:);
    TrainData1=Data(1:trainno,:);
    XValidation1=Data(trainno+1:trainno+valno,:);
    Responses1=categorical(ones(trainno,1),cl);
    YValidation1=categorical(ones(valno,1),cl);

    xx=xx+1;
    Data=DataREST2;
    shuffling=randperm(task_nsubn{xx});
    Data=Data(shuffling(1:task_nsubn{xx}),:);
    TrainData2=Data(1:trainno,:);
    XValidation2=Data(trainno+1:trainno+valno,:);
    Responses2=categorical(2*ones(trainno,1),cl);
    YValidation2=categorical(2*ones(valno,1),cl);

    TrainData=vertcat(TrainData1,TrainData2);
    XValidation=vertcat(XValidation1,XValidation2);
    Responses=vertcat(Responses1,Responses2);
    YValidation=vertcat(YValidation1,YValidation2);
        
    t = templateSVM('KernelFunction','polynomial');
    svmmodel=fitcecoc(TrainData,Responses,'Learners',t);

    %% compute 
    
    con=zeros(2,2);
    test1=predict(svmmodel,XValidation1);
    for i=1:valno
        winclass=test1(i);
        con(1,winclass)=con(1,winclass)+1;
    end
    test2=predict(svmmodel,XValidation2);
    for i=1:valno
        winclass=test2(i);
        con(2,winclass)=con(2,winclass)+1;
    end
    
    con=con/valno;
    accdist(nfold,:)=diag(con);
    pc=pc+con;
end
pc=pc/kfold
acc=sum(diag(pc))/2

save results_classMOVIE_1vs4.mat acc pc accdist;
