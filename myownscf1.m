clear all
close all
fc=200;%载波频率 
fs=2000;%采样频率
fd=80;%波特率
N=fs/fd;%每个符号的采样点数
M=800;%符号数
ts=1/fs;%采样时间间隔
t0=10.5;
t=0:ts:t0;
N1=200;
N2=M*N/N1;%总采样点数/N1
select=menu('调制方式','2PSK','4PSK','2FSK','4FSK','MSK');
switch select
    case 1,%2PSK调制
        x=randint(1,M);%随机产生元素为0,1长度为M的向量，概率各为1/2
        y=zeros(1,M*N);
        for i=1:M
                for j=1:N
                    y((i-1)*N+j)=cos(2*pi*fc*t((i-1)*N+j)+(1-x(i))*pi);
                end
        end
        %y=wgn(1,M*N,0);
    case 2,%4PSK调制
        x=randint(1,M,4);%随机产生元素为0,1，2,3长度为M的向量，概率各为1/4
        y=zeros(1,M*N);
        for i=1:M
                for j=1:N
                    y((i-1)*N+j)=cos(2*pi*fc*t((i-1)*N+j)+pi/2*x(i));
                end
        end
    case 3,%2FSK调制
        x=randint(1,M);  
        y=zeros(1,M*N);
        for i=1:M
               for j=1:N;
                    f0=fc-2*fd+4*x(i)*fd;
                    y((i-1)*N+j)=cos(2*pi*f0*t((i-1)*N+j));
                end
           
        end
    case 4,%4FSK调制
        x=randint(1,M,4);
        y=zeros(1,M*N);
        for i=1:M
               for j=1:N
                   if x(i)<2
                       f0=fc-(x(i)+1)*fd;
                       y((i-1)*N+j)=cos(2*pi*f0*t((i-1)*N+j));
                   else
                       f0=fc+(x(i)-1)*fd;
                       y((i-1)*N+j)=cos(2*pi*f0*t((i-1)*N+j));
                   end
               end
        end
    case 5,%MSK调制
         D=randint(1,M);
        d1=2*D-1;
        D1=[];
        x=ones(1,N);
        D1=[];
        for n=1:M
            D1=[D1,d1(n)*x];
        end   
        fai=zeros(1,M*N);
        fai(1:N)=0;
        for l=2:M
            for i=1:N
                if d1(l)==d1(l-1)
                    fai((l-1)*N+1:l*N)=fai((l-2)*N+1:(l-1)*N);
                else
                    fai((l-1)*N+1:l*N)=fai((l-2)*N+1:(l-1)*N)+(l-1)*pi;
                end
            end
        end
        for i=1:N*M
            y(i)=cos(2*pi*(fc+D1(i)*fd/4)*t(i)+fai(i));
        end
end
    
signal2=zeros(N2,N1);

for  j=1:N2
     for i=1:N1
           signal2(j,i)=y(1,N1*(j-1)+i);
     end
end

 spectrum=zeros(N2,N1);
for i=1:N2
 spectrum(i,:)=fft(signal2(i,:));%进行离散傅里叶变换
 spectrum(i,:)=fftshift(spectrum(i,:));%将前半部分元素与后半部分元素交换
end

scf=zeros(N1,N1/2);
alfa=-N1/2:2:N1/2-2;
%f=0:N1-1;
for  i=1:N2
   for f=0:N1-1
        for h=0:N1/2-1
            d=mod(f+alfa(h+1)/2,N1);%Circular shift
            d=d+1;%The index of a matrix starts at 1, but d is from 0 to N1-1, so d=d+1
            d2=mod(f-alfa(h+1)/2,N1);%Circular shift
            d2=d2+1;
            scf(f+1,h+1)=scf(f+1,h+1)+spectrum(i,d2)*spectrum(i,d)';
        end
   end
end
n=(-N1/2:2:N1/2-2)/N1;
m=(-N1/2:1:N1/2-1)/N1;
%m=(N1-1:-1:0)/N1
scf1=abs(scf);
p=max(max(scf1));
[X,Y]=meshgrid(n,m);%确定X,Y范围
figure(1)
mesh(X,Y,scf1/p);
%title('循环自相关归一化函数图');
xlabel('α / Fsamp');
ylabel('f / Fsamp');
zlabel('magnitude');
scf=max(scf1);
figure(2)
stem(n,scf/p)