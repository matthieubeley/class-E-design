%This script computes first the normalized class E circuit parameters.
%Then it computes the actual circuit parameters, considering the operating point
%(Switching frequency f, Input DC voltage V and RF Load resistance R)

%The input are :
%   q --> the input LC tank resonant frequency ratio
%   D --> the ON duty cycle
%   K --> the switch voltage at turn-on (0 if ZVS)
%   s --> the switch voltage derivative at turn-on (0 is ZdVS, or ZCS)

%Please use at the maximum two of these four inputs as vectors using
%a = linspace(a_min,a_max,number of points a_i)
%and the other transposed using
%b = linspace(b_min,b_max,number of points b_i).'
%and the other inputs as scalars


%% Normalized circuit parameters

% Design point (§§§ TO ADAPT §§§)
D = 0.5;            % Duty cycle
% D = linspace(0.1,0.9,50).';
q = 1.2915;         % Frequency ratio (!= 0 and != 1 becuase of numerical issues)
% q = linspace(0.01,2.7,500);
K = 0;              % Normalized switch voltage at turn-on
% K = linspace(0,1,11).';
s = 0;              % Normalized switch voltage derivative at turn-on
% s = linspace(-2,2,50);


% Analytical expressions for the normalized circuit (DO NOT MODIFY)
a1 = q./(q.^2-1).*(cos(q*2*pi.*(D-1))-q.^2.*cos(2*pi*D)+(q.^2-1));
b1 = q./(q.^2-1).*(-q.*sin(q*2*pi.*(D-1))+(q.^2).*sin(2*pi*D));
c1 = 2*pi.*q.*D+(K-1).*sin(2*pi*q.*(D-1))+s./q.*(1-cos(2*pi*q.*(D-1)));
a2 = q./(q.^2-1).*(sin(q*2*pi.*(D-1))-q.*sin(2*pi*D));
b2 = q./(q.^2-1).*(q.*cos(q*2*pi.*(D-1))-q.*cos(2*pi*D));
c2 = -(K-1).*cos(2*q*pi.*(D-1))-1-s./q.*sin(2*pi*q.*(D-1));

paramA = (b1.*c2-b2.*c1)./(a2.*b1-a1.*b2);
paramB = (a1.*c2-a2.*c1)./(a2.*b1-a1.*b2);
p = sqrt(paramA.^2+paramB.^2);
phi = atan2(paramA,paramB);

k1 = p.*q./(1-q.^2).*(q.*cos(2*pi*q).*cos(phi)+sin(2*pi*q).*sin(phi))+(K-1).*cos(2*pi*q)-s./q.*sin(2*pi*q);
k2 = p.*q./(1-q.^2).*(q.*sin(2*pi*q).*cos(phi)-cos(2*pi*q).*sin(phi))+(K-1).*sin(2*pi*q)+s./q.*cos(2*pi*q);

g = (pi*D.^2)./p-cos(2*pi*D+phi)/(2*pi)-D.*sin(phi)+cos(phi)/(2*pi)+D./q.^2./p.*s;

X_q = 1./(q.^2-1).*(k1.*(cos(2*pi*q.*D).*sin(2*pi*D+phi)-cos(2*pi*q).*sin(phi)-q.*sin(2*pi*q.*D).*cos(2*pi*D+phi)+q.*sin(2*pi*q).*cos(phi)) ...
    +k2.*(sin(2*pi*q.*D).*sin(2*pi*D+phi)-sin(2*q*pi).*sin(phi)+q.*cos(2*pi*q.*D).*cos(2*pi*D+phi)-q.*cos(2*pi*q).*cos(phi)));
X_0 = sin(phi)-sin(2*pi*D+phi);
X_n = q.^2./(q.^2-1).*p/4.*(4*pi*(1-D)+sin(2*phi)-sin(4*pi*D+2*phi));
R_q = 1./(q.^2-1).*(k1.*(-cos(2*pi*q.*D).*cos(2*pi*D+phi)+cos(2*pi*q).*cos(phi)-q.*sin(2*pi*q.*D).*sin(2*pi*D+phi)+q.*sin(2*pi*q).*sin(phi)) ...
    +k2.*(-sin(2*pi*q.*D).*cos(2*pi*D+phi)+sin(2*pi*q).*cos(phi)+q.*cos(2*pi*q.*D).*sin(2*pi*D+phi)-q.*cos(2*pi*q).*sin(phi)));
R_0 = -cos(phi)+cos(2*pi*D+phi);
R_n = q.^2./(q.^2-1).*p/4.*(-cos(2*phi)+cos(4*pi*D+2*phi));

% Switch voltage first harmonic
Vs1X_norm = 1/pi*(X_q+X_0+X_n);             
Vs1R_norm = 1/pi*(R_q+R_0+R_n);            
Vs1_norm = sqrt(Vs1X_norm.^2+Vs1R_norm.^2); 
psi = atan(Vs1X_norm./Vs1R_norm);      


% Normalized circuit parameters (DO NOT MODIFY)
L_norm = p/2./g; % Input inductance
C_norm = 2*g./p./q.^2; % Shunt capacitance
X_norm = -Vs1X_norm./Vs1R_norm; % Additional series reactance
P_norm = 2*g.^2; % Power
Idc_norm = 2*g.^2; % Input DC current
Vout_norm = 2*g; % Output voltage magnitude
Rdc_norm = 1/2./g.^2; % Input equivalent resistor
Ir_norm = 2*g; % Output current magnitude


%Switch constraints

Vsp_norm = (1.7613+0.05*q)./(1-D); % Peak switch voltage

Isp_norm = 2*g.*(2*pi*D./p - sin(phi) + sin(2*pi*D+phi)); % Peak switch current
for Di = 1:length(D)
    for qi = 1:length(q)
        if ((p(Di,qi) >= 1) && (2*pi*D(Di)/p(Di,qi)+sin(2*pi*D(Di)+phi(Di,qi))-sin(acos(-1/p(Di,qi)))-(acos(-1/p(Di,qi))-phi(Di,qi))/p(Di,qi) <= 0) && (acos(-1/p(Di,qi))-phi(Di,qi)-2*pi*D(Di) <= 0))
            Isp_norm(Di,qi) = 2*g(Di,qi)*((acos(-1/p(Di,qi))-phi(Di,qi))/p(Di,qi)-sin(phi(Di,qi))+sin(acos(-1/p(Di,qi))));
        end
    end
end

% RMS switch current
IsRMS_norm = 2*g./p.*sqrt(8/6*pi^2*D.^3-2*pi*D.^2.*p.*sin(phi)+p.^2/pi.*sin(phi).*(cos(2*pi*D+phi)-cos(phi))-2*D.*p.*cos(2*pi*D+phi)-p.^2.*sin(4*pi*D+2*phi)./8/pi+p/pi.*sin(2*pi*D+phi)-D.*p.^2.*cos(2*phi)/2+D.*p.^2-p/pi.*sin(phi)+p.^2.*sin(2*phi)/8/pi);

% Power output capability P/(Isp*Vsp)
cp = 2*g.^2./Vsp_norm./Isp_norm;

% Modified power output capability P/(IsRMS*Vsp)
cpmr = 2*g.^2./Vsp_norm./IsRMS_norm;


%% Denormalized circuit parameters

% Operating point (§§§ TO ADAPT §§§)
f = 6.78e6;     % Switching frequency [Hz]
V = 30;         % Input DC voltage [V]
R = 20;         % Load resistance [Ohm]

% Derived parameters
w = 2*pi*f;     % Switching pulsation [rad/s]

% Circuit parameters
L = L_norm*R/w;
C = C_norm/R/w;
X = X_norm*R;
Rdc = Rdc_norm*R;
P = P_norm*V^2/R;
Idc = Idc_norm*V/R;
Vout = Vout_norm*V;
Ir = Ir_norm*V/R;
Vs1 = Vs1_norm*V;
Vs1R = Vs1R_norm*V;
Vs1X = Vs1X_norm*V;
Vsp = Vsp_norm*V;
Isp = Isp_norm*V/R;
IsRMS = IsRMS_norm*V/R;


% Filter design
Q = 7;         % Quality factor, should be high enough (§§§ TO ADAPT §§§)
L0 = Q*R/w;     % Filter inductance
C0 = 1/(Q*R*w); % Filter capacitance
% Series reactance X should be considered (added to L0 or C0)



%% Theoretical normalized waveforms (Vin = 1 V, R = 1 Ohm, w = 1 rad/s) 
% Only with scalar inputs q, D, K and s
% Comment out if q, D, K, and s are not all scalars

% Time vector (one period)
T = 2*pi;
n_points = 5000;
Tsample = T/n_points;
t = (0:n_points-1)*Tsample;

% Switch voltage
vs = k1*cos(q*t)+k2*sin(q*t)+1+q^2/(q^2-1)*p*cos(t+phi);
vs(t<(D*T))=0;

% Switch voltage first harmonic
vs1 = Vs1R_norm*sin(t+phi)+Vs1X_norm*cos(t+phi);

% Output voltage
vout = -2*g*sin(t+phi);

% Output resonant branch current
ir = 2*g*sin(t+phi);

% Input current
iin = 2*g/p*(t-p*sin(phi)+s/q^2);
iin(t>(D*T)) = 2*g/p/q*(-k1*sin(q*t(t>(D*T)))+k2*cos(q*t(t>(D*T)))-p*q^3/(q^2-1)*sin(t(t>(D*T))+phi));

% Switch current
is = iin+ir;
is(t>(D*T)) = 0;

% Shunt capacitance current
ic = iin+ir; 
ic(t<(D*T))=0;