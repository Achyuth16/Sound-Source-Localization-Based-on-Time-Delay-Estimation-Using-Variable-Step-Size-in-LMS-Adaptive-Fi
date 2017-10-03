clc;
clear all;
close all;
x1=load('data1.mat');
x2=load('data2.mat');
x3=load('data3.mat');
stem3(x1.xes,x1.yes,x1.zes,'r');
hold on;grid on;
stem3(x1.mx,x1.my,x1.mz);
hold on
stem3(x2.xes,x2.yes,x2.zes,'r');
hold on;
stem3(x2.mx,x2.my,x2.mz);
hold on;
stem3(x3.xes,x3.yes,x3.zes,'r');
hold on;
stem3(x3.mx,x3.my,x3.mz);
xlabel('xlabel/m');
ylabel('ylabel/m');
zlabel('zlabel/m');
legend('estimate','theoretical')
% axis([-50,50,-50,50,-50,50]);
e1=sqrt((x3.mx-x3.xes)^2+(x3.my-x3.yes)^2+(x3.mz-x3.zes)^2);