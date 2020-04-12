%Athanasiou Ioannis
%First Task

clc;
clear all;

SNRdb=1:1:20; %SNR values in dB from 1 to 20
SNR=10.^(SNRdb/10);
mon_car=10^4; %number of loops monte carlo
Eb=1;
mesh_timh_ch=1; %mean E[h]=1
diaspora_ch=1; %variance of h=   1
T=10; % gia ka8e T bits to h paramenei sta8ero
num_bits=input('num_bits= ');
transmission_rate=input('transmission_rate= ');

No=2*Eb./SNR; %Noise variance
diasp=No/2;

errorsV = zeros(1,length(diasp));
errorsY = errorsV;
N=zeros(length(diasp), num_bits);
Y=N;
V=Y;



for k=1:1:mon_car
    s_bits=randi([0 1],1,num_bits); %random stream with zeros and ones
    h=diaspora_ch.*randn(1,floor(num_bits/T)+1)+mesh_timh_ch;%random variable h
    
    %oversampling
    os_h=upsample(h,T);
	for i=1:T:length(os_h)
        for j=1:T-1
            os_h(i+j)=os_h(i);
        end
    end
    
	os_h=os_h(1:num_bits);
    stream=2*s_bits-1; %convert to stream with 1 and -1
    
    for z=1:1:length(diasp)
        N(z,:)=sqrt(diasp(z))*randn(1, num_bits); 
        Y(z,:)=os_h.*stream+N(z,:);
        %(EQUALIZATION)
		V(z,:)=os_h.^(-1).*Y(z,:);
        for m=1:num_bits
            %check for correct or wrong guess
            if ((V(z,m)<0) && (stream(m)==1))||((V(z,m)>0) && (stream(m)==-1))
                errorsV(z)=errorsV(z)+1;
            end
            if ((Y(z,m)<0) && (stream(m)==1))||((Y(z,m)>0) && (stream(m)==-1))
                errorsY(z)=errorsY(z)+1;
            end
        end
    end
    %scatterplot(Y);
end


err_loopV=errorsV/mon_car;
BER_simV=err_loopV/num_bits;

err_loopY=errorsY/mon_car;
BER_simY=err_loopY/num_bits;

BER_theor=qfunc(sqrt(SNR));
%graph SNRvsBER
semilogy(SNRdb,BER_simV,'k');
hold on
semilogy(SNRdb,BER_simY,'b-.');
semilogy(SNRdb,BER_theor,'r*');
legend('Simulation with equalization','Simulation without equalization','Theoritical',3);
xlabel('SNR in dB');
ylabel('BER');
axis([min(SNRdb) max(SNRdb) 10^(-5) 1]);
title('BER vs SNR');
hold off


disp('Press Enter gia to ERWTHMA A.b');

pause;


%Diversity


Div_br=input('branches= ');

errors= zeros(1,length(diasp));
N=zeros(length(diasp), num_bits);
Y=N;
V=Y;
h=zeros(Div_br,floor(num_bits/T)+1);


for k=1:1:mon_car
    os_h=zeros(Div_br,length(h)*T);
    s_bits=randi([0 1],1,num_bits);%random stream with zeros and ones
    stream=2*s_bits-1; %convert to stream with 1 and -1
    for i_ant=1:Div_br
        h(i_ant,:)=diaspora_ch.*randn(1,floor(num_bits/T)+1)+mesh_timh_ch;
        os_h(i_ant,:)=upsample(h(i_ant,:),T);
    end
	for i_ant=1:Div_br
        for i=1:T:length(os_h(i_ant,:))
            for j=1:T-1
                os_h(i_ant,i+j)=os_h(i_ant,i);
            end
        end
	end
    [r,c]=size(os_h);
	diaf=c-num_bits;

    while diaf>0
		os_h(:,c)=[];
		[r,c]=size(os_h);
		diaf=diaf-1;
    end
	os_h_div=os_h.^2;
    par=sum(os_h_div,1);
    for z=1:1:length(diasp)
		for i_ant=1:Div_br
			N(z,:,i_ant)=sqrt(diasp(z))*randn(1, num_bits); 
			Y(z,:,i_ant)=os_h(i_ant,:).*stream+N(z,:,i_ant);
            V(z,:,i_ant)=os_h(i_ant,:).*Y(z,:,i_ant);
        end
        % MRC according to formula (h1y1+h2y2...)/(h1^2+h2^2....)
        U=sum(V,3);
        MRC=U.*(par(1).^(-1));
        for m=1:num_bits
            %check for correct or wrong guess
            if ((MRC(z,m)<0) && (stream(m)==1))||((MRC(z,m)>0) && (stream(m)==-1))
                errors(z)=errors(z)+1;
            end
        end
    end
end

err_loop=errors/mon_car;
BER_sim=err_loop/num_bits;


BER_theor=qfunc(sqrt(SNR));
%graph SNRvsBER
semilogy(SNRdb,BER_sim,'k');
hold on
semilogy(SNRdb,BER_theor,'r*');
legend('Simulation with MRC','Theoritical',3);
xlabel('SNR in dB');
ylabel('BER');
axis([min(SNRdb) max(SNRdb) 10^(-5) 1]);
title('BER vs SNR with MRC');
hold off

disp('Press Enter for ERWTHMA B')
pause;



%Block Coding

errorsV = zeros(1,length(diasp));
errorsY = errorsV;
N=zeros(length(diasp), num_bits);
Y=N;
V=Y;
h=zeros(1,Div_br);


for k=1:1:mon_car
    s_bits=randi([0 1],1,num_bits); %random stream with zeros and ones
    h=diaspora_ch.*randn(1,floor(num_bits/T)+1)+mesh_timh_ch;
    os_h=upsample(h,T);
	for i=1:T:length(os_h)
        for j=1:T-1
            os_h(i+j)=os_h(i);
        end
    end
	os_h=os_h(1:num_bits);
    stream=2*s_bits-1; %convert to stream with 1 and -1
    
    for z=1:1:length(diasp)
        N(z,:)=sqrt(diasp(z))*randn(1, num_bits); 
        Y(z,:)=os_h.*stream+N(z,:); 
		V(z,:)=(os_h+1).^(-1).*Y(z,:);
        for m=1:num_bits
            %check for correct or wrong guess
            if ((V(z,m)<0) && (stream(m)==1))||((V(z,m)>0) && (stream(m)==-1))
                errorsV(z)=errorsV(z)+1;
            end
            if ((Y(z,m)<0) && (stream(m)==1))||((Y(z,m)>0) && (stream(m)==-1))
                errorsY(z)=errorsY(z)+1;
            end
        end
    end
    %scatterplot(Y);
end

err_loopV=errorsV/mon_car;
BER_simV=err_loopV/num_bits;
err_loopY=errorsY/mon_car;
BER_simY=err_loopY/num_bits;
n = 23;             % Codeword length
k = 12;             % Message length
dmin = 7;           % Minimum distance
codedBER = bercoding(SNRdb,'block','hard',n,k,dmin); 

BER_theor=qfunc(sqrt(SNR));
%graph SNRvsBER
semilogy(SNRdb,BER_simV,'k');
hold on
semilogy(SNRdb,BER_simY,'b-.');
semilogy(SNRdb,BER_theor,'r*');
semilogy(SNRdb,codedBER,'m--');
legend('Simulation with equalization','Simulation without equalization','Theoritical','Coded BER',3);
xlabel('SNR in dB');
ylabel('BER');
axis([min(SNRdb) max(SNRdb) 10^(-5) 1]);
title('BER vs SNR for h=h+1');
hold off
pause;

%Second Task

clc;
clear all;

SNRdb=1:1:20; %SNR values in dB from 1 to 20
SNR=10.^(SNRdb/10);
mon_car=10^4; %number of loops monte carlo
Eb=1;
mesh_timh_ch=1; %mean E[h]=1
diaspora_ch=1; %variance of h=   1
T=10; % 
num_bits=input('num_bits= ');
transmission_rate=input('transmission_rate= ');

No=2*Eb./SNR; %Noise variance
diasp=No/2;


in_ant=2;
out_ant=2;


errors= zeros(1,length(diasp));
N=zeros(length(diasp), num_bits);
Y=N;
V=Y;
h=zeros(Div_br,floor(num_bits/T)+1);


for k=1:1:mon_car
    os_h=zeros(Div_br,length(h)*T);
    s_bits=randi([0 1],1,num_bits); %random stream with zeros and ones
    stream=2*s_bits-1; %convert to stream with 1 and -1
    for i_ant=1:in_ant+out_ant
        h(i_ant,:)=diaspora_ch.*randn(1,floor(num_bits/T)+1)+mesh_timh_ch;
        os_h(i_ant,:)=upsample(h(i_ant,:),T);
    end
	for i_ant=1:in_ant+out_ant
        for i=1:T:length(os_h(i_ant,:))
            for j=1:T-1
                os_h(i_ant,i+j)=os_h(i_ant,i);
            end
        end
	end
    [r,c]=size(os_h);
	diaf=c-num_bits;

    while diaf>0
		os_h(:,c)=[];
		[r,c]=size(os_h);
		diaf=diaf-1;
    end
    for z=1:1:length(diasp)
		for i_ant=1:in_ant
			N(z,:,i_ant)=sqrt(diasp(z))*randn(1, num_bits); 
			Y(z,:,i_ant)=os_h(i_ant,:).*stream+N(z,:,i_ant); 
        end
        for m=1:num_bits
            %check for correct or wrong guess
            if ((MRC(z,m)<0) && (stream(m)==1))||((MRC(z,m)>0) && (stream(m)==-1))
                errors(z)=errors(z)+1;
            end
        end
    end
end

err_loop=errors/mon_car;
BER_sim=err_loop/num_bits;


BER_theor=qfunc(sqrt(SNR));
%graph SNRvsBER
semilogy(SNRdb,BER_sim,'k');
hold on
semilogy(SNRdb,BER_theor,'r*');
legend('Simulation with spatial multiplexing for a system 2x2','Theoritical',3);
xlabel('SNR in dB');
ylabel('BER');
axis([min(SNRdb) max(SNRdb) 10^(-5) 1]);
title('BER vs SNR with spatial multiplexing for a system 2x2');
hold off

pause;

%Third Task

clc;
clear all;

SNRdb=1:1:20; %SNR values in dB from 1 to 20
SNR=10.^(SNRdb/10); 
mon_car=10^4; %number of loops monte carlo
Eb=1;
mesh_timh_ch=1; %mean E[h]=1
diaspora_ch=1; %variance of h=   1
T=10;
num_bits=250;
transmission_rate=1;

No=2*Eb./SNR; %Noise variance
diasp=No/2;


errorsV = zeros(1,length(diasp)); %size of errors for question B
errorsY = errorsV; %size of errors for question B
N=zeros(num_bits+4,length(diasp)); %size of noise for question B
Y=N; %size of output y (from y=hx+n)
V=zeros(num_bits,length(diasp)); %size for frequency equalizer
H=zeros(num_bits+4,num_bits); %size of array H.




for k=1:1:mon_car
    counter=1;
    s_bits=randi([0 1],1,num_bits); %random stream with zeros and ones
    h=diaspora_ch.*randn(1,5)+mesh_timh_ch;%random variable h
    for i=1:1:num_bits+4
        for j=1:1:num_bits
            if mod(counter,5)==0
                h=diaspora_ch.*randn(1,5)+mesh_timh_ch;
            end
            if i==j
                H(i,j)=h(1,1);
            elseif i-j==1
                H(i,j)=h(1,2);
            elseif i-j==2
                H(i,j)=h(1,3);
            elseif i-j==3
                H(i,j)=h(1,4);
            elseif i-j==4
                H(i,j)=h(1,5);
            end
            counter=counter+1;
        end
    end
            
    
    stream=2*s_bits-1; %convert to stream with 1 and -1
    
    for z=1:1:length(diasp)
        N(:,z)=sqrt(diasp(z))*randn(num_bits+4,1);
        Y(:,z)=H*stream'+N(:,z); 
        %(EQUALIZATION)
		V(:,z)=(((H'*H)^(-1))*H')*Y(:,z);
        for m=1:num_bits
            %check for correct or wrong guess
            if ((V(m,z)<0) && (stream(m)==1))||((V(m,z)>0) && (stream(m)==-1))
                errorsV(z)=errorsV(z)+1;
            end
            if ((Y(m,z)<0) && (stream(m)==1))||((Y(m,z)>0) && (stream(m)==-1))
                errorsY(z)=errorsY(z)+1;
            end
        end
    end
end

err_loopV=errorsV/mon_car;
BER_simV=err_loopV/num_bits;

err_loopY=errorsY/mon_car;
BER_simY=err_loopY/num_bits;

%BER_theor=qfunc(sqrt(SNR));
%graph SNRvsBER
semilogy(SNRdb,BER_simV,'k');
hold on
semilogy(SNRdb,BER_simY,'b-.');
%semilogy(SNRdb,BER_theor,'r*');
legend('Simulation for Question B','Simulation for Question A',3);
xlabel('SNR in dB');
ylabel('BER');
axis([min(SNRdb) max(SNRdb) 10^(-5) 1]);
title('BER vs SNR');
hold off
