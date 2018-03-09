% % % Skript pro výpočet změny U, H, S, F a G pro 1 mol chlóru

% % Nadefinované podmínky
T1 = 298; % K
p1 = 01*10^6; % Pa

T2 = 500; % K
p2 = 10*10^6; % Pa

% % Tabelované hodnoty a funkce (webbook.nist.gov/cgi/cbook.cgi?ID=C7782505&Mask=1)
% (použitelné pouze v rozmezí 298-1000 K)
% Konstanty:
A = 33.05060;
B = 12.22940;
C = -12.06510;
D = 4.385330;
E = -0.159494;
G = 259.0290;

% Funkce
% St = @(t) A*log(t) + B*t + (C*t^2)/2 + (D*t^3)/3 -E/(2*t^2) + G;
Cpot = @(t) A + B*t + C*t.^2 + D*t.^3 - E./(t.^2);
% S = @(t) Sot(t/1000); % J/(mol*K)
Cpo = @(t) Cpot(t/1000); % J/(mol*K)

% Hodnoty
So = 223.08; % J/(mol*K)
S1 = So;

% % Mezivýpočty
% Přepočet Cp na Cv
R = 8314;
Cvo = @(t) Cpo(t) - R;

% Výpočet V1 a V2
V1 = R*T1/p1;
V2 = R*T2/p2;

% % Výpočet delS
funkce_ = @(t) Cpo(t)./t;
S2 = S1 + 0 + integral(funkce_, T1, T2) - R*log(p2/p1) + 0; % 0 znamenají nulové integrály R/p* - R/p, protože počítáme plyn ve stavu IP
delS = S2 - S1;

% % Výpočet delU
funkce_ = @(t) Cvo(t);
delU = 0 + integral(funkce_, T1, T2) + 0;

% % Výpočet delH
funkce_ = @(t) Cpo(t);
delH = 0 + integral(funkce_, T1, T2) + 0;

% % Výpočet delG
delG = delH - T2*S2 + T1*S1;

% % Výpočet delF
delF = delU - T2*S2 + T1*S1;

% Výsledné výpočty
delS
fprintf('   J/(mol*K)\n')
delU
fprintf('   J\n')
delH
fprintf('   J\n')
delG
fprintf('   J\n')
delF
fprintf('   J\n')
