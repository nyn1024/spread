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
    spectrum(i,:)=fft(signal2(i,:));%进行离散傅里叶变换
    spectrum(i,:)=fftshift(spectrum(i,:));%将前半部分元素与后半部分元素交换
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
cycret=max_a/max(max_a);   % 得到每个alfa对应的最大值(归一化)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% 
n=(-N1/2:2:N1/2-2)/N1;     %坐标归一化
m=(-N1/2:1:N1/2-1)/N1;
% 
% 
% %size(scf1)    %??????
p=max(max(scf));
% [X,Y]=meshgrid(n,m);%确定X,Y范围,          网格采样点的函数
% figure(1)
% mesh(X,Y,scf/p);  %  绘制由线条框构成的曲面
% title('循环自相关归一化函数图');
% xlabel('循环频率alfa（fs)');
% ylabel('频率f (fs）')
scf=max(scf);
% figure(2)
stem(n,scf/p)    %stem函数用于绘制火柴梗图

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


