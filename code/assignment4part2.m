%% Assignemnt 4 Part 2: Transient Circuit Simulation
%% Part 2: Time Domain MNA Techniques

% a) This is a low pass filter. It remains a linear circuit. 
% B) We would expect low frequency signals to be transmitted with a gain,
% and higher freqeuncy signals to be attenuated.

%% Parameters
clc
G1=1/1;
C2=0.25;
G2=0.5;
L=0.2;
G3=0.1;
G4=100;
G5=1/1000;
ALPHA=100;
Vin=10;
vx=-10;

%% Constructing the C, G, and F matrices

%X=[V1 Iin V2 V3 V4 V5 IL I4]
G=[-G1,1, G1, 0,0,0,0,0; ...
    G1,0, -G1-G2,0,0,0,-1,0;...
    0,0,0,-G3,0,0,1,0;...
    0,0,1,-1,0,0,0,0;...
    0,0,0,0,0,G4,-ALPHA*G4,1;...
    0,0,0,0,0,-G4-G5,ALPHA*G4,0;...
    1,0,0,0,0,0,0,0;
    0,0,0,0,1,0,-ALPHA,0];

C = [-C2,0,C2,0,0,0,0,0;...
    C2,0,-C2,0,0,0,0,0;....
    0,0,0,0,0,0,0,0;...
    0,0,0,0,0,0,-L,0;...
    0,0,0,0,0,0,0,0;...
    0,0,0,0,0,0,0,0;...
    0,0,0,0,0,0,0,0;...
    0,0,0,0,0,0,0,0];

%% MATLAB simulation 

%time step
deltaT=0.001;

A=C/(deltaT)+G;

%% 1. Unit Step Simulation
 Vjm=[0; 0; 0; 0; 0; 0; 0; 0];

voutval=zeros(1000,1);
vinval2=zeros(1000,1);
count=1;
time=zeros(1000,1);


for t=0:deltaT:1
    if(t>0.03)
        F=[0 0 0 0 0 0 1 0];
    else
        F=[0 0 0 0 0 0 0 0];
    end
    
    time(count)=t;
    Vj=inv(A)*(C*Vjm/deltaT + F');
    voutval(count)=Vj(6);
    vinval2(count)=Vj(1);
    Vjm=Vj;
    count=count+1; 
end

figure(1)
subplot(2,1,2)
plot(time,voutval)
title('Unit Step: Vout over Time')
grid on

subplot(2,1,1)
plot(time,vinval2)
title('Unit Step: Vin over Time')
grid on

%Vout Frequency plots

fs=1000;
fvout=fft(voutval);
n=length(voutval);

Y=fftshift(fvout);
fshift = (-n/2:n/2-1)*(fs/n); % zero-centered frequency range
powershift = abs(Y).^2/n;     % zero-centered power

figure(2)
plot(fshift,powershift)
title('Vout: Unit Step F Spectrum')


%Vin

fvin=fft(vinval2);
n=length(vinval2);

Y=fftshift(fvin);
fshift = (-n/2:n/2-1)*(fs/n); % zero-centered frequency range
powershift = abs(Y).^2/n;     % zero-centered power

figure(3)
plot(fshift,powershift)
title('Vin: Unit Step F Spectrum')

%% 2. Sine Function

vx=@(t) sin(2*pi*33 *t);
Vjm=[0; 0; 0; 0; 0; 0; 0; 0];
voutval=zeros(1000,1);
vinval2=zeros(1000,1);
count=1;
time=zeros(1000,1);
    
    
     for t=0:deltaT:1
    
         time(count)=t;
         F=[0 0 0 0 0 0 vx(t) 0];
         Vj=inv(A)*(C*Vjm/deltaT + F');
         voutval(count)=Vj(6);
         vinval2(count)=Vj(1);
         Vjm=Vj;
         count=count+1;
    
      
     end
     
     figure(4)
     subplot(2,1,2)
     plot(time,voutval)
     title('Sine Function: Vout over Time')
     grid on
     
     subplot(2,1,1)
     plot(time,vinval2)
     title('Sine Function: Vin over Time')
     grid on
     
%Vout Frequency plots
  
fs=1000;
fvout=fft(voutval);
n=length(voutval);
Y=fftshift(fvout);
fshift = (-n/2:n/2-1)*(fs/n); % zero-centered frequency range
powershift = abs(Y).^2 /n;     % zero-centered power

figure(5)
plot(fshift,powershift)
title('Vout: Sine F Spectrum')

%Vin

fvin=fft(vinval2);
n=length(vinval2);
Y=fftshift(fvin);
fshift = (-n/2:n/2-1)*(fs/n); % zero-centered frequency range
powershift = abs(Y).^2/n;     % zero-centered power

figure(6)
plot(fshift,powershift)
title('Vin: Sine F Spectrum')

%% 3. Gaussian Function

     vx=@(t) exp(-(1/2)*((t-0.06)/(0.03))^2);
     Vjm=[0; 0; 0; 0; 0; 0; 0; 0];
     voutval=zeros(1000,1);
     vinval2=zeros(1000,1);
     count=1;
     time=zeros(1000,1);
    
     
     for t=0:deltaT:1
    
         time(count)=t;
         F=[0 0 0 0 0 0 vx(t) 0];
         Vj=inv(A)*(C*Vjm/deltaT + F');
         voutval(count)=Vj(6);
         vinval2(count)=Vj(1);
         Vjm=Vj;
         count=count+1;
         
        
     end
     
     figure(7)
    
     subplot(2,1,2)
     plot(time,voutval)
     title('Guassian Function: Vout over Time')
     grid on
     
     subplot(2,1,1)
     plot(time,vinval2)
     title('Gaussian Function: Vin over Time')
     grid on
     
%Vout Frequency plots
fs=1000;
fvout=fft(voutval);
n=length(voutval);
Y=fftshift(fvout);
 fshift = (-n/2:n/2-1)*(fs/n); % zero-centered frequency range
 powershift = abs(Y).^2/n;     % zero-centered power

 figure(8)
plot(fshift,powershift)
title('Vout: Gaussian F Spectrum')


%Vin

fvin=fft(vinval2);
n=length(vinval2);
t=(0:n-1)*(fs/n);
power = abs(fvin).^2 / n;
 Y=fftshift(fvin);
 fshift = (-n/2:n/2-1)*(fs/n); % zero-centered frequency range
 powershift = abs(Y).^2/n;     % zero-centered power

 figure(9)
plot(fshift,powershift)
title('Vin: Gaussian F Spectrum')      
        

%% 4. Sine Function (Small and Large Time Steps)

% smaller time step first
deltaT=0.001;

A=C/(deltaT)+G;
vx=@(t) sin(2*pi*33 *t);
Vjm=[0; 0; 0; 0; 0; 0; 0; 0];
voutvalsmall=zeros(1000,1);
vinval2small=zeros(1000,1);
count=1;
timesmall=zeros(1000,1);
    
    
     for t=0:deltaT:1
    
         timesmall(count)=t;
         F=[0 0 0 0 0 0 vx(t) 0];
         Vj=inv(A)*(C*Vjm/deltaT + F');
         voutvalsmall(count)=Vj(6);
         vinval2small(count)=Vj(1);
         Vjm=Vj;
         count=count+1;
    
      
     end
 
% larger time step

deltaT=0.01;

A=C/(deltaT)+G;

vx=@(t) sin(2*pi*33 *t);
Vjm=[0; 0; 0; 0; 0; 0; 0; 0];
voutvallarge=zeros(100,1);
vinval2large=zeros(100,1);
count=1;
timelarge=zeros(100,1);
    
    
     for t=0:deltaT:1
    
         timelarge(count)=t;
         F=[0 0 0 0 0 0 vx(t) 0];
         Vj=inv(A)*(C*Vjm/deltaT + F');
         voutvallarge(count)=Vj(6);
         vinval2large(count)=Vj(1);
         Vjm=Vj;
         count=count+1;
    
      
     end
  figure(10)
     subplot(2,1,2)
     plot(timesmall,voutvalsmall)
     title('Sine Function: Vout over Small Time')
     grid on
     
     subplot(2,1,1)
     plot(timelarge,voutvallarge)
     title('Sine Function: Vin over Large Time')
     grid on
%% 5. Conclusion
% The time domain response of the circuit was observed in this module of
% the report. It can be observed in part 4 that the larger the time step
% the more distorted the output waveform appears. Therefore a smaller time
% step should be used for this analysis. 