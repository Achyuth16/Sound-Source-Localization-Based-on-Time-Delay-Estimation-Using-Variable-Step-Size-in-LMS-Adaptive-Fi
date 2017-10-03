clc;
clear all;
close all;
%%
% din=randn(1,1000);
R=50;% distance between source and reference microphone
varphi=120/180*pi;% angle of direction(-180--180)
theta=120/180*pi;% angle of dip(0--180)
D=5;% distance between microphones
C=343.39;% speed of sound

mx=R*sin(theta)*cos(varphi); % theoretical x,y,z
my=R*sin(theta)*sin(varphi);
mz=R*cos(theta);

mr1=sqrt(mx^2+my^2+mz^2);%  mathmatical model 
md10=sqrt((mx-D)^2+my^2+mz^2)-mr1;
md20=sqrt(mx^2+(my-D)^2+mz^2)-mr1;
md30=sqrt((mx+D)^2+my^2+mz^2)-mr1;
md40=sqrt(mx^2+(my+D)^2+mz^2)-mr1;
md50=sqrt(mx^2+my^2+(mz-D)^2)-mr1;
% md60=sqrt(mx^2+my^2+(mz+D)^2)-mr1;

mt10=md10/C;% convert to theoretical time delay
mt20=md20/C;
mt30=md30/C;
mt40=md40/C;
mt50=md50/C;
% mt60=mt60/C;

fs=10000; % sampling frequency
t=0:1/fs:4-1/fs;% N=40000;
din=chirp(t,500,1,2000);
% din=randn(1,1000);
md=round(fs*[mt10 mt20 mt30 mt40 mt50]);% transfer to delay of samples
dmax=max(md);
x0=[zeros(1,2*dmax-1),din(1:length(din)-(2*dmax-1))];% theoretical received signals without noise
x1=[zeros(1,2*dmax-md(1)-1),din(1:length(din)-(2*dmax-md(1)-1))];
x2=[zeros(1,2*dmax-md(2)-1),din(1:length(din)-(2*dmax-md(2)-1))];
x3=[zeros(1,2*dmax-md(3)-1),din(1:length(din)-(2*dmax-md(3)-1))];
x4=[zeros(1,2*dmax-md(4)-1),din(1:length(din)-(2*dmax-md(4)-1))];
x5=[zeros(1,2*dmax-md(5)-1),din(1:length(din)-(2*dmax-md(5)-1))];

snr=5;%SNR
% powdin=sum(din.^2)/length(din)*fs;
powdin=sum(din.^2)/length(din);
pownoise=powdin/(10^(snr/10));

x0=x0+sqrt(pownoise)*randn(1,length(din));% with noise
x1=x1+sqrt(pownoise)*randn(1,length(din));
x2=x2+sqrt(pownoise)*randn(1,length(din));
x3=x3+sqrt(pownoise)*randn(1,length(din));
x4=x4+sqrt(pownoise)*randn(1,length(din));
x5=x5+sqrt(pownoise)*randn(1,length(din));

N=length(x0);
k=320;% filter order
b=[zeros(1,k/2),1];%delay L/2 of the desired signal
a=1;
x1=filter(b,a,x1);
x2=filter(b,a,x2);
x3=filter(b,a,x3);
x4=filter(b,a,x4);
x5=filter(b,a,x5);


for a=1:5
    if (a==1)
        dn=x1;
    elseif (a==2)
            dn=x2;
    elseif (a==3)
                dn=x3;
    elseif (a==4)
                    dn=x4;
    elseif(a==5)
                        dn=x5;
    end  
    
%% LMS 
% q=0.0005; % step size 
% h=zeros(1,k)'; 
% for i=k:N  
%     
%     u=x0(i:-1:i-k+1); 
%     en(i-k+1)=dn(i)-u*h;   
%     h=h+q*en(i-k+1)*u';%update weight
% 
% end;

%% variable step size LMS
beta=0.00114;%Variable step size factor
alpha=0.00114;%Variable step size factor
h=zeros(1,k)';%update vector
for i=k:N  
    
    u=x0(i:-1:i-k+1);
    en(i-k+1)=dn(i)-u*h;   
    q=beta*(1-exp(-alpha*(en(i-k+1)^2)));%update step size
    h=h+q*en(i-k+1)*u';%update weight

end;
 
% stem(h);
num=h;
den=1;
wf=freqz(num,den);
phasewf=phase(wf);
printph(:,a)=phase(wf);% phase
f=(0:length(wf)-1)*fs/length(wf);
ki=(phasewf(500)-phasewf(200))/(pi*(f(500)-f(200))/fs);%calculate derivation k
DT(a)=ki+k/2;%find the delay in samples
Tn0(a)=DT(a)/fs;%the dalays in sec

end
plot(f,printph(:,1)*180/pi);
hold on;grid on;
plot(f,printph(:,2)*180/pi,'r');
hold on;
plot(f,printph(:,3)*180/pi,'k');
hold on;
plot(f,printph(:,4)*180/pi,'g');
hold on;
plot(f,printph(:,5)*180/pi,'--m');
xlabel('Frequency(Hz)')
ylabel('Phase(degree)')
legend('t10','t20','t30','t40','t50');

t10=Tn0(1);
t20=Tn0(2);
t30=Tn0(3);
t40=Tn0(4);
t50=Tn0(5);

r=(4*D^2-C^2*(t10^2+t20^2+t30^2+t40^2))/(2*(t10+t20+t30+t40)*C);% get r
gvarphi=atan(((t40-t20)*(t40*C+t20*C+2*r))/((t30-t10)*(t30*C+t10*C+2*r)))/pi*180;% get varphi
gtheta=asin(sqrt(((t40*C-t20*C)*(t40*C+t20*C+2*r))^2+((t30*C-t10*C)*(t30*C+t10*C+2*r))^2)/(r*4*D))/pi*180;%get theta

if(t20<0&&t30<0)
    gvarphi=180+gvarphi;
elseif(t30<0&&t40<0)
    gvarphi=-180+gvarphi;
end

if (t50>0)% determine theta is larger than 90 or not
    gtheta=180-gtheta;
end
xes=r*sin(gtheta/180*pi)*cos(gvarphi/180*pi);%calculate eastimated x
yes=r*sin(gtheta/180*pi)*sin(gvarphi/180*pi);%calculate eastimated y
zes=r*cos(gtheta/180*pi);% calculate eastimated z

e1=abs(mt10-t10)/abs(mt10)*100;
e2=abs(mt20-t20)/abs(mt20)*100;
e3=abs(mt30-t30)/abs(mt30)*100;
e4=abs(mt40-t40)/abs(mt40)*100;
e5=abs(mt50-t50)/abs(mt50)*100;
% error=sqrt((xes-mx)^2+(yes-my)^2+(zes-mz)^2)
figure(2)
stem3(xes,yes,zes);
hold on;grid on;
stem3(mx,my,mz,'r');
axis([-50,50,-50,50,-50,50]);
save data3 mx my mz xes yes zes







































