%NiDUC-Przesylanie informacji z wykorzystaniem modulacji PSK
nS = 100;
nSym = 1024; %liczba bitow generacji
%szum = 30
M = 2;
Tb = 1e-6;
fc = 1e6; %Czestotliwosc nosna
T=1/fc;
t2=T/99:T/99:T;
 
%GENERATOR BITOW******************************************************
s = randi([0 M-1], nSym, 1);
%s_noise = awgn(s,20);
 
szum = input('Prosze podac wspolczynnik szumu: ');
%BPSK*****************************************************************
s_bsk = pskmod(s, M, pi);
rx_bsk = awgn(s_bsk,20,szum);
s_bsk = rectpulse(s_bsk, nS);
%figure('Name', 'Kanal BPSK');
h1 = scatterplot(rx_bsk);
hold on
title('BPSK');
rx_bsk2 = rectpulse(rx_bsk, nS);
bit=[[0 M-1], nSym, 1];
for i=1:length(rx_bsk)
    if(rx_bsk(i)>0)
        bit(i)=1;
    elseif(rx_bsk(i)<0)
        bit(i)=0;
    end
end
w = pskmod(bit, M, pi);
bit_BSK = bit.';
%w = pskmod(bit, M, pi);
w2 = rectpulse(w, nS);

%QPSK*****************************************************************
x = -1.5:.01:1.5;
xz = -x;
data_NZR=2*s-1;
s_p_data=reshape(data_NZR,2,length(s)/2);
y=[];
y_in=[];
y_qd=[];
d=zeros(1,length(s)/2);
for i=1:length(s)/2
    p=s(2*i);
    imp=s(2*i - 1);
    y1=s_p_data(1,i)*cos(2*pi*fc*t2); % inphase component
    y2=s_p_data(2,i)*sin(2*pi*fc*t2) ;% Quadrature component
    y_in=[y_in y1]; % inphase signal vector
    y_qd=[y_qd y2]; %quadrature signal vector
    y=[y (y1+y2)*0.7]; % modulated signal vector
    if (imp == 0) && (p == 0)
       d(i)=exp(pi/4);%45 degrees
    end
    if (imp == 1)&&(p == 0)
        d(i)=exp(3*pi/4);%135 degrees
    end
    if (imp == 1)&&(p == 1)
        d(i)=exp(5*pi/4);%225 degrees
    end
    if (imp == 0)&&(p == 1)
        d(i)=exp(7*pi/4);%315 degrees
    end
end
Tx_sig=y; % transmitting signal after modulation
qpsk=d;
tt=T/99:T/99:(T*length(s))/2;
 
%KANAL***************************************************************
s_re = reshape(s,2,nSym/2);
s_2 = zeros(nSym/2, 1);
for i=1:nSym/2
    if(s_re(1,i) == 0 && s_re(2,i) == 0)
        s_2(i,1) = 1;
    elseif(s_re(1,i) == 0 && s_re(2,i) == 1)
        s_2(i,1) = 2;    
    elseif(s_re(1,i) == 1 && s_re(2,i) == 0)
        s_2(i,1) = 0;        
    else
        s_2(i,1) = 3;
    end      
end
s_qpsk = pskmod(s_2, 2*M, pi); %Non-Return to zero level encoder
rx_qpsk = awgn(s_qpsk,20,szum);
s_qpsk2 = rectpulse(s_qpsk, nS);
%figure('Name', 'Kanal QPSK');
h2 = scatterplot(rx_qpsk);
hold on
title('QPSK');

bit2=[[0 M-1], nSym, 1];
bit2=[[0 M-1], nSym, 1];
for i=1:length(rx_qpsk)
    if((angle(rx_qpsk(i))>2.3562 && angle(rx_qpsk(i))<3.1416) || (angle(rx_qpsk(i))<-2.3562 && angle(rx_qpsk(i))>-3.1416))
        bit2(2*i-1)=1;
        bit2(2*i)=0;
    elseif(angle(rx_qpsk(i))>0.7854 && angle(rx_qpsk(i))<2.3562)
        bit2(2*i-1)=1;
        bit2(2*i)=1;
    elseif(angle(rx_qpsk(i))<-0.7854 && angle(rx_qpsk(i))>-2.3562)
        bit2(2*i-1)=0;
        bit2(2*i)=0;
    elseif((angle(rx_qpsk(i))>-0.7854 && angle(rx_qpsk(i))<0) || (angle(rx_qpsk(i))>0 && angle(rx_qpsk(i))<0.7854))
        bit2(2*i-1)=0;
        bit2(2*i)=1;
    end
end

wq = pskmod(bit2, M, pi);
%w = pskmod(bit, M, pi);
bit_QPSK = bit2.';
wq2 = rectpulse(wq, nS);
 
%CZAS****************************************************************
t = 0: (Tb / nS) :nSym*Tb - (Tb / nS); %Time Domain
t = transpose(t);
 
%WIZUALIZACJA WYNIKOW************************************************
figure('Name', 'Kanal');
 
subplot(1,1,1);
plot(t, rx_bsk2,'m','linewidth',2);
title('Sygnal z szumem');
xlabel('Czas');
ylabel('Amplituda');
axis([0 (nSym*Tb) -1.2 1.2]);
 
figure('Name', 'Generacja sygnalu i Modulacja BPSK i QPSK')
 
subplot(3,3,1);
plot(t, s_bsk,'linewidth',2);
axis([0 (nSym*Tb) -1.2 1.2]); %Plot ONLY first 10 bits
title('Bity wejsciowe');
xlabel('Czas');
ylabel('Amplituda');
 
s_c_bsk = s_bsk.*cos(2*pi*fc*t);
subplot(3,1,2);
plot(t, s_c_bsk,'r', 'linewidth',2);
axis([0 (nSym*Tb) -1.2 1.2]);
title('Modulacja BSK');
xlabel('Czas');
ylabel('Amplituda');
 
subplot(3,1,3);
plot(tt,Tx_sig,'g','linewidth',2), grid on;
title('Modulacja QPSK');
xlabel('Czas');
ylabel('Amplituda');
axis([0 (nSym*Tb) -1.2 1.2]);
 
subplot(3,3,2);
plot(t, w2,'linewidth',2);
axis([0 (nSym*Tb) -1.2 1.2]); %Plot ONLY first 10 bits
title('Demodulacja BSK');
xlabel('Czas');
ylabel('Amplituda');

subplot(3,3,3);
plot(t, wq2,'linewidth',2);
axis([0 (nSym*Tb) -1.2 1.2]); %Plot ONLY first 10 bits
title('Demodulacja QPSK');
xlabel('Czas');
ylabel('Amplituda');