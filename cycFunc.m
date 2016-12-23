function cycret=cycFunc(signal_input)

N1=200;
N2=length(signal_input)/N1;
signal2=zeros(N2,N1);


for  j=1:N2
    for i=1:N1
        signal2(j,i)=signal_input(1,N1*(j-1)+i);
    end
end

spectrum=zeros(N2,N1);
for i=1:N2
    spectrum(i,:)=fft(signal2(i,:));%������ɢ����Ҷ�任
    spectrum(i,:)=fftshift(spectrum(i,:));%��ǰ�벿��Ԫ�����벿��Ԫ�ؽ���
end

scf=zeros(N1,N1/2);
alfa=-N1/2:2:N1/2-2;

for  i=1:N2    %1~100
    for f=0:N1-1             %0~199        
        for h=0:N1/2-1      %0~99
            d=mod(f+alfa(h+1)/2,N1);%Circular shift
            d=d+1;%The index of a matrix starts at 1, but d is from 0 to N1-1, so d=d+1
            d2=mod(f-alfa(h+1)/2,N1);%Circular shift
            d2=d2+1;
            scf(f+1,h+1)=scf(f+1,h+1)+spectrum(i,d2)*spectrum(i,d)';
        end
    end
end
scf=abs(scf);
[max_a,index]=max(scf);
cycret=max_a/max(max_a);   % �õ�ÿ��alfa��Ӧ�����ֵ(��һ��)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% 
n=(-N1/2:2:N1/2-2)/N1;     %�����һ��
m=(-N1/2:1:N1/2-1)/N1;
% 
% 
% %size(scf1)    %??????
p=max(max(scf));
% [X,Y]=meshgrid(n,m);%ȷ��X,Y��Χ,          ���������ĺ���
% figure(1)
% mesh(X,Y,scf/p);  %  �����������򹹳ɵ�����
% title('ѭ������ع�һ������ͼ');
% xlabel('ѭ��Ƶ��alfa��fs)');
% ylabel('Ƶ��f (fs��')
scf=max(scf);
% figure(2)
stem(n,scf/p)    %stem�������ڻ��ƻ��ͼ

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


