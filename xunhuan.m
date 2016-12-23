clc;
close all;

fs=5000;%码率
Ts=1/fs;%码片周期
fc=1000;%载频
threshold=10e-16;

SNR=0;
number_of_code=1e3;
modulation_type=menu('调制方式','QPSK','16QAM','64QAM');
switch modulation_type
    case 1;%QPSK
        load s_PSK;
      signal=s_PSK; 
%  signal_send=Windows_send(s_PSK,100);%信号过窗(2PSK,4PSk需要过窗)*****这里需要注意，2PSK和4PSK需要，其他***
%  signal=transpose(signal_send);
%  t=(Ts)*[1:length(signal)];
%  y=real(signal).*cos(2*pi*fc*t)-imag(signal).*sin(2*pi*fc*t);
%  cyc_ret=cycFunc(y(1:2000));
%  x = -0.5:1/100:0.5-0.01;
%  stem(x,cyc_ret);
    case 2;%16QAM
        load s_16QAM;
        signal=s_16QAM;
        
    case 3;%64QAM
        load s_64QAM;
        signal=s_64QAM;
end     
       
%QPSK信号生成        
%  M=4;
%  number_of_code=1e4;
%  phase_offset=0;
%  code=randi([0 M-1],number_of_code,1);
%  PSK_mod=comm.PSKModulator('ModulationOrder',M);
%  PSK_Demod=comm.PSKDemodulator('ModulationOrder',M);
%  PSK_mod.PhaseOffset = phase_offset;
%  s_PSK=step(PSK_mod,code);
% save s_PSK.mat;
%  title_str=[num2str(M),'PSK'];
%  scatter(real(signal),imag(signal));
%QAM信号生成
% M=16;
% number_of_code=1e4;
% code=randi([0 M-1],number_of_code,1);
% QAM_mod=comm.RectangularQAMModulator('ModulationOrder',M);
% QAM_Demod=comm.RectangularQAMDemodulator('ModulationOrder',M);
% QAM_mod.PhaseOffset = 0;
% s_16QAM=step(QAM_mod,code);
% save s_16QAM.mat;

%%%%%扩频前信号的循环谱图%%%%%%%%%%%%%%%%%
%  t=(Ts)*[1:length(signal)];
%  y=real(signal).*cos(2*pi*fc*t)-imag(signal).*sin(2*pi*fc*t);
%  cyc_ret=cycFunc(y(1:2000));
%  x = -0.5:1/100:0.5-0.01;
%  stem(x,cyc_ret); 
% 
%%%%%%扩频和解扩过程%%%%%%%%%%%%%%%% 
load spreading_code.mat;
spreaded_code_matrix=kron(signal,spreading_code);%%%%%%%%%%%扩频%%%%%%%%%%%%signal是n行1列
% [m,n]=size(temp);
% signal_spreading=reshape(transpose(temp),1,m*n);
% % signal_send=Windows_send(signal_spreading1,100);%信号过窗(2PSK,4PSk需要过窗)*****这里需要注意，2PSK和4PSK需要，其他****
% signal_spreading=transpose(signal_send);
% signal_spreading_length=length(signal_spreading);
% y=real(signal_spreading).*cos(2*pi*fc*t)-imag(signal_spreading).*sin(2*pi*fc*t);
% signal_spreading=awgn(signal_spreading,SNR,'measured');
despread_length=4;
despread_number=1;
despreading_code=ovsf(despread_length);%%%%产生解扩码%%%%
despreading_code=despreading_code(despread_number,:);
% signal_spreading2=Windows_receive(signal_spreading,100);
% temp1=reshape(signal_spreading,numel(signal_spreading)/despread_length,despread_length);
% spreaded_code_matrix=transpose(temp1);
spreaded_code_sequence=reshape(transpose(spreaded_code_matrix),numel(spreaded_code_matrix),1);
spreaded_code_matrix2=reshape(spreaded_code_sequence,despread_length,numel(spreaded_code_matrix)./despread_length);
spreaded_code_matrix2=transpose(spreaded_code_matrix2);
despreading_code_matrix=repmat(despreading_code,size(spreaded_code_matrix2,1),1);
signal_despreading=(despreading_code_matrix.*spreaded_code_matrix2)./despread_length;
signal_despreading=sum(signal_despreading,2);
signal_despreading=(signal_despreading>=threshold).*signal_despreading;
% signal_despreading=transpose(signal_despreading);
signal_send1=Windows_send(signal_despreading,100);%信号过窗(2PSK,4PSk需要过窗)*****这里需要注意，2PSK和4PSK需要，其他***
signal_despreading1=transpose(signal_send1);
 t=(Ts)*[1:length(signal_despreading1)];
 y=real(signal_despreading1).*cos(2*pi*fc*t)-imag(signal_despreading1).*sin(2*pi*fc*t);
 cyc_ret=cycFunc(y(1:2000));
 x = -0.5:1/100:0.5-0.01;
 stem(x,cyc_ret); 
xlabel('循环频率');
title('解扩码16-1');