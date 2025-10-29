%% Inputs

%Duty cycle
%D = round(linspace(0.1,0.9,351),5).'; %Column vector
D = 0.5; %Scalar

%q ratio
%q = round(linspace(0,4,800),5); %Row vector
q = 1.412;   %q != 1   %Scalar

%Operating point
Vin = 1; %Input DC voltage
f = 1/2/pi; %Switching frequency
RL = 1; %Load resistance

w = 2*pi*f; %Angular frequency
T = 1/f; %Period

%% Design of the normalized class E inverter (RL = 1 Ohm, w = 1 rad/s, Vin = 1 V)

%Constant parameters calculation
a1 = q./(q.^2-1).*(cos(q*2*pi.*(D-1))-q.^2.*cos(2*pi*D)+q.^2-1);
b1 = q./(q.^2-1).*(-q.*sin(q*2*pi.*(D-1))+q.^2.*sin(2*pi*D));
c1 = 2*pi.*q.*D-sin(2*pi*q.*(D-1));
a2 = q./(q.^2-1).*(sin(q*2*pi.*(D-1))-q.*sin(2*pi*D));
b2 = q./(q.^2-1).*(q.*cos(q*2*pi.*(D-1))-q.*cos(2*pi*D));
c2 = cos(2*q*pi.*(D-1))-1;

paramA = (b1.*c2-b2.*c1)./(a2.*b1-a1.*b2);
paramB = (a1.*c2-a2.*c1)./(a2.*b1-a1.*b2);
p = sqrt(paramA.^2+paramB.^2);
phi = atan2(paramA,paramB);

k1 = p.*q./(1-q.^2).*(q.*cos(2*pi*q).*cos(phi)+sin(2*pi*q).*sin(phi))-cos(2*pi*q);
k2 = p.*q./(1-q.^2).*(q.*sin(2*pi*q).*cos(phi)-cos(2*pi*q).*sin(phi))-sin(2*pi*q);

gx = (pi*D.^2)./p-cos(2*pi*D+phi)/(2*pi)-D.*sin(phi)+cos(phi)/(2*pi);

%Analytical derivation of the parameters for the switch voltage First Harmonic Approximation (FHA)
X1 = 1./(q.^2-1).*(k1.*(cos(2*pi*q.*D).*sin(2*pi*D+phi)-cos(2*pi*q).*sin(phi)-q.*sin(2*pi*q.*D).*cos(2*pi*D+phi)+q.*sin(2*pi*q).*cos(phi)) ...
    +k2.*(sin(2*pi*q.*D).*sin(2*pi*D+phi)-sin(2*q*pi).*sin(phi)+q.*cos(2*pi*q.*D).*cos(2*pi*D+phi)-q.*cos(2*pi*q).*cos(phi)));
X2 = sin(phi)-sin(2*pi*D+phi);
X3 = q.^2./(q.^2-1).*p/4.*(4*pi*(1-D)+sin(2*phi)-sin(4*pi*D+2*phi));
R1 = 1./(q.^2-1).*(k1.*(-cos(2*pi*q.*D).*cos(2*pi*D+phi)+cos(2*pi*q).*cos(phi)-q.*sin(2*pi*q.*D).*sin(2*pi*D+phi)+q.*sin(2*pi*q).*sin(phi)) ...
    +k2.*(-sin(2*pi*q.*D).*cos(2*pi*D+phi)+sin(2*pi*q).*cos(phi)+q.*cos(2*pi*q.*D).*sin(2*pi*D+phi)-q.*cos(2*pi*q).*sin(phi)));
R2 = -cos(phi)+cos(2*pi*D+phi);
R3 = q.^2./(q.^2-1).*p/4.*(-cos(2*phi)+cos(4*pi*D+2*phi));

VX_norm = 1/pi*(X1+X2+X3);
VR_norm = 1/pi*(R1+R2+R3);
psi = atan(VX_norm./VR_norm);

%Analytical derivation of the parameters for the input current spectral decomposition (up to harmonic #n_harm)
n_harm = 10; %Number of harmonics
harm_i = zeros(1,1,n_harm);
harm_i(1,1,:) = linspace(1,n_harm,n_harm);

Mu1 = 2*pi*D./harm_i.*sin(2*pi*D.*harm_i)+(cos(2*pi*D.*harm_i)-1)./(harm_i.^2)+q-q;
Mu2 = -p./harm_i.*sin(phi).*sin(2*pi*D.*harm_i);
Mu3 = -k1./q./(harm_i.^2-q.^2).*(-sin(2*pi*D.*q).*sin(2*pi*D.*harm_i).*harm_i+q.*cos(2*pi.*q)-q.*cos(2*pi*D.*q).*cos(2*pi*D.*harm_i));
Mu4 = k2./q./(harm_i.^2-q.^2).*(-cos(2*pi*D.*q).*sin(2*pi*D.*harm_i).*harm_i-q.*sin(2*pi.*q)+q.*sin(2*pi*D.*q).*cos(2*pi*D.*harm_i));
Mu5 = -p.*q.^2./((q.^2-1).*(harm_i.^2-1)).*(cos(phi)-harm_i.*sin(2*pi*D+phi).*sin(2*pi*D.*harm_i)-cos(2*pi*D+phi).*cos(2*pi*D.*harm_i));
Mu5(:,:,1) = -p.*q.^2./(q.^2-1).*(-cos(phi)/4+pi*(1-D).*sin(phi)+cos(4*pi*D+phi)/4);
Mun = 2*gx./(pi*p).*(Mu1+Mu2+Mu3+Mu4+Mu5);

Nu1 = -2*pi*D./harm_i.*cos(2*pi*D.*harm_i)+sin(2*pi*D.*harm_i)./(harm_i.^2)+q-q;
Nu2 = p.*sin(phi).*(cos(2*pi*D.*harm_i)-1)./harm_i;
Nu3 = -k1./q./(harm_i.^2-q.^2).*(-harm_i.*sin(2*pi*q)+harm_i.*sin(2*pi*D.*q).*cos(2*pi*D.*harm_i)-q.*cos(2*pi*q.*D).*sin(2*pi*D.*harm_i));
Nu4 = k2./q./(harm_i.^2-q.^2).*(-harm_i.*cos(2*pi*q)+harm_i.*cos(2*pi*D.*q).*cos(2*pi*D.*harm_i)+q.*sin(2*pi*q.*D).*sin(2*pi*D.*harm_i));
Nu5 = -p.*q.^2./((q.^2-1).*(harm_i.^2-1)).*(-harm_i.*sin(phi)+harm_i.*sin(2*pi*D+phi).*cos(2*pi*D.*harm_i)-cos(2*pi*D+phi).*sin(2*pi*D.*harm_i));
Nu5(:,:,1) = -p.*q.^2./(q.^2-1).*(-sin(phi)/4+pi*(1-D).*cos(phi)+sin(4*pi*D+phi)/4);
Nun = 2*gx./(pi*p).*(Nu1+Nu2+Nu3+Nu4+Nu5);

%Analytical derivation of the parameters for the input current RMS value
A1 = 4/3*pi^2*D.^3;
A2 = -2*pi.*D.^2.*p.*sin(phi);
A3 = D.*p.^2.*sin(phi).*sin(phi);
B1 = 1/2/pi*(k1./q).^2.*(pi*(1-D)+1/4./q.*(sin(4*pi*q.*D)-sin(4*pi*q)));
B2 = 1/2/pi*(k2./q).^2.*(pi*(1-D)+1/4./q.*(sin(4*pi*q)-sin(4*pi*q.*D)));
B3 = 1/2/pi*(p.*q.^2./(q.^2-1)).^2.*(pi*(1-D)+1/4*(sin(4*pi*D+2*phi)-sin(2*phi)));
B4 = k1.*k2/4/pi./(q.^3).*(cos(4*pi*q)-cos(4*pi*q.*D));
B5 = k1.*p.*q/pi./((q.^2-1).^2).*(-sin(2*pi*q.*D).*cos(2*pi*D+phi)+sin(2*pi*q).*cos(phi)-q.*cos(2*pi*q).*sin(phi)+q.*cos(2*pi*q.*D).*sin(2*pi*D+phi));
B6 = k2.*p.*q/pi./((q.^2-1).^2).*(cos(2*pi*q.*D).*cos(2*pi*D+phi)-cos(2*pi*q).*cos(phi)-q.*sin(2*pi*q).*sin(phi)+q.*sin(2*pi*q.*D).*sin(2*pi*D+phi));
paramC = sqrt(A1+A2+A3+B1+B2+B3+B4+B5+B6);

%Normalized circuit parameters
Lp_norm = p/2./gx;
Cp_norm = 2*gx./p./q.^2;
X_norm = VX_norm./VR_norm;
P_norm = 2*gx.^2;

%Normalized voltage and currents
Iin_norm = 2*gx.^2;
Vout_norm = 2*gx;
Rdc_norm = 1/2./gx.^2;
Ir_norm = 2*gx;
Vs_fund_norm = sqrt(VX_norm.^2+VR_norm.^2);
Iin_RMS_norm = 2*gx./p.*paramC;

%Spectral decomposition (3rd dimension) of input current
Iin_spectr_norm = zeros(length(D),length(q),n_harm+1);
Iin_spectr_norm(:,:,1) = Iin_norm;
Iin_spectr_norm(:,:,2:end)=sqrt(Mun.^2+Nun.^2);

%Normalized maximum switch voltage
sV = (1.7613+0.05*q)./(1-D);

%Normalized maximum switch current
sI = 2*gx.*(2*pi*D./p - sin(phi) + sin(2*pi*D+phi));
for Di = 1:length(D)
    for qi = 1:length(q)
        if ((p(Di,qi) >= 1) && (2*pi*D(Di)/p(Di,qi)+sin(2*pi*D(Di)+phi(Di,qi))-sin(acos(-1/p(Di,qi)))-(acos(-1/p(Di,qi))-phi(Di,qi))/p(Di,qi) <= 0) && (acos(-1/p(Di,qi))-phi(Di,qi)-2*pi*D(Di) <= 0))
            sI(Di,qi) = 2*gx(Di,qi)*((acos(-1/p(Di,qi))-phi(Di,qi))/p(Di,qi)-sin(phi(Di,qi))+sin(acos(-1/p(Di,qi))));
        end
    end
end

%Normalized RMS switch current
sIRMS = 2*gx./p.*sqrt(8/6*pi^2*D.^3-2*pi*D.^2.*p.*sin(phi)+p.^2/pi.*sin(phi).*(cos(2*pi*D+phi)-cos(phi))-2*D.*p.*cos(2*pi*D+phi)-p.^2.*sin(4*pi*D+2*phi)./8/pi+p/pi.*sin(2*pi*D+phi)-D.*p.^2.*cos(2*phi)/2+D.*p.^2-p/pi.*sin(phi)+p.^2.*sin(2*phi)/8/pi);

%Power output capability
cp = 2*gx.^2./sV./sI; %with Imax
cpmr = 2*gx.^2./sV./sIRMS; %with Irms
 
%% Design of the de-normalized class E inverter

%Circuit parameters
Lp = Lp_norm * RL/w;
Cp = Cp_norm * 1/RL/w;
X = X_norm * RL;
P = P_norm * Vin^2/RL;
Rdc = Rdc_norm * RL;

%Circuit voltages and currents
Iin = Iin_norm * Vin/RL;
Vout = Vout_norm * Vin;
Ir = Ir_norm * Vin/RL;
Vs_fund = Vs_fund_norm * Vin;
Iin_RMS = Iin_RMS_norm * Vin/RL;
Iin_spectr = Iin_spectr_norm * Vin/RL;

%Switch constraints
Vp = sV * Vin;
Ip = sI * Vin/RL;
IsRMS = sIRMS * Vin/RL;

%% Theoretical waveforms (only once q and D are scalar !!)

%Time vector
n_per = 3; %Period number
npoints_per = 1000; %Number of points per period
fsample = npoints_per*f; %Sample frequency
Tsample = 1/fsample; %Sample period
npoints = n_per*npoints_per; %Number of points
t = (0:npoints-1)*Tsample; %Time vector
t1 = t(1:npoints_per); %Time vector (one period)

%Drive voltage vgs
Vdrive = 5; %Vdd drive voltage
vgs1 = zeros(1,npoints_per);
vgs1(t1<(D*T)) = Vdrive; %One period
vgs = repmat(vgs1,1,n_per); %n_per periods

%Switch voltage vs
vs1 = Vin * (k1*cos(q*t1*w)+k2*sin(q*t1*w)+1+q^2/(q^2-1)*p*cos(t1*w+phi));
vs1(t1<(D*T))=0; %One period
vs = repmat(vs1,1,n_per); %n_per periods

%Fundamental of switch voltage vs1
vs_fund1 = Vin * (VR_norm*sin(t1*w+phi)+VX_norm*cos(t1*w+phi)); %One period
vs_fund = repmat(vs_fund1,1,n_per); %n_per periods

%Output voltage vout
vout1 = -Vin * 2*gx*sin(t1*w+phi); %One period
vout = repmat(vout1,1,n_per); %n_per periods

%Output current ir
ir1 = Vin/RL * 2*gx*sin(t1*w+phi); %One period
ir = repmat(ir1,1,n_per); %n_per periods

%Input current iin
iin1 = Vin/RL * 2*gx/p*(t1*w-p*sin(phi));
iin1(t1>(D*T)) = Vin/RL * 2*gx/p*(-k1/q*sin(q*t1(t1>(D*T))*w)+k2/q*cos(q*t1(t1>(D*T))*w)-p*q^2/(q^2-1)*sin(t1(t1>(D*T))*w+phi)); %One period
iin = repmat(iin1,1,n_per); %n_per periods

%Switch current is
is1 = iin1+ir1;
is1(t1>(D*T)) = 0; %One period
is = repmat(is1,1,n_per); %n_per periods

%Shunt capacitor current ic
ic1 = iin1+ir1; 
ic1(t1<=(D*T))=0; %One period
ic = repmat(ic1,1,n_per); %n_per periods
