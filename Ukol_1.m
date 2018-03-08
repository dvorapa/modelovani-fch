% Skript pro výpočet změny U, H, S, F a G pro 1 mol chlóru

% Podmínky
T1 = 298; % K
p1 = 0.1*10^6; % Pa

T2 = 500; % K
p2 = 10*10^6; % Pa

% Data a funkce z tabulek (webbook.nist.gov/cgi/cbook.cgi?ID=C7782505&Mask=1)
% použitelné pouze v rozmezí 298-1000 K
A = 33.05060;
B = 12.22940;
C = -12.06510;
D = 4.385330;
E = -0.159494;
G = 259.0290;
So = 223.08; % J/(mol*K)
Sot = @(t) A*log(t) + B*t + (C*t^2)/2 + (D*t^3)/3 -E/(2*t^2) + G;
Cpot = @(t) A + B*t + C*t^2 + D*t^3 - E/(t^2);
So = @(t) Sot(t/1000); % J/(mol*K)
Cpo = @(t) Cpot(t/1000); % J/(mol*K)

% Přepočet Cp na Cv
R = 8.314;
Cvo = @(t) Cpo(t) - R;

% vlastní výpočty
delU = 0;
delH = 0;
delS = 0;
delF = 0;
delG = 0;