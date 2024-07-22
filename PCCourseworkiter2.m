close all
clear all
clc
% Введение переменных
Na = 6.022e26; %с учетом перевода в кг
M_CO = 28.01; 
M_H2O = 18.01528; 
M_CO2 = 44.01; 
M_H2 = 2.016;
M_C = 12.011; 
M_H = 1.00794; 
M_O = 15.999;
M = [28.01;18.01528;44.01;2.016];
% g/mol
m_CO = M_CO/Na; 
m_H2O = M_H2O/Na; 
m_CO2 = M_CO2/Na; 
m_H2 = M_H/Na; 
m_C = M_C/Na; 
m_H = M_H/Na; 
m_O = M_O/Na;
% kg
r_H2 = 0.74e-10; 
r_CO = 1.128e-10; 
r_CO2 = 1.16e-10; 
r_H2O = 9.57e-11;
% m
angle_H2O = 1.8242;
% rad
T = [298 500 1000];
% K
R = 8.31446261815324;
% J/(mol*K)
k = 1.380649e-23; 
P = 101325;
% Pa
c = 299792458;
% m*s^2
h = 6.62607015e-34;
% kg*m^2/sec
I_CO = (m_C*m_O)/(m_C+m_O)*r_CO^2; 
I1_H2O = (2*m_O*m_H)/(2*m_H+m_O)*r_H2O^2*(cos(angle_H2O/2))^2; 
I2_H2O = 2*m_H*r_H2O^2*(sin(angle_H2O/2))^2; 
I3_H2O = I1_H2O+I2_H2O; 
I_H2O = [I1_H2O; I2_H2O; I3_H2O]; 
I_CO2 = (2*m_O)*r_CO2^2; 
I_H2 = (1/2)*m_H*r_H2^2;
I = [I_CO; I_H2O; I_CO2; I_H2];
% значения вычисленны при помощи стр. 29-31 в методических рекомендациях.
sigma_CO = 1; 
sigma_H2O = 1; 
sigma_CO2 = 2; 
sigma_H2 = 2;
sigma = [1;1;2;2];
% коэффициент симметрии с учетом геометрии молекулы.
theta_CO = h*c*2169e2/k; 
theta1_H2O = h*c*3656e2/k; 
theta2_H2O = h*c*3755.8e2/k; 
theta3_H2O = h*c*1594.8e2/k; 
theta_H2O = [theta1_H2O;theta2_H2O;theta3_H2O]; 
theta1_CO2 = h*c*1388.17e2/k; 
theta2_CO2 = h*c*2349.16e2/k; 
theta3_CO2 = h*c*667.4e2/k; 
theta_CO2 = [theta1_CO2; theta2_CO2;theta3_CO2]; 
theta_H2 = h*c*4396.554e2/k;
theta = [theta_CO; theta_H2O; theta_CO2; theta_H2];
% характеристические частоты найдены при помощи стр. 183-189 КСФХВ, так же
% от туда были взяты весовые коэффициенты вырождения колебаний.
G = [1;1;1;1];
% все молекулы имеют терм Сигма

% расчет поступательной составляющей.
lnQtrans = (1.5).*log(M)+(2.5).*log(T)-log(P)+8.8612; 
datalnQtrans = array2table(lnQtrans,'VariableNames',{'298K' '500K' '1000K'}, 'RowNames',{'CO', 'H2O', 'CO2','H2'});
Utrans = 1.5*R*T; 
Htrans = 2.5*R*T;
HH0trans = 2.5*R*T;
Gibbstrans = 20.786.*log(T)+12.4716.*log(M)-30.4760; 
Strans = 12.4716.*log(M)+20.7860.*log(T)-9.6853; 
dataStrans = array2table(Strans,'VariableNames',{'298K' '500K' '1000K'}, 'RowNames',{'CO', 'H2O', 'CO2','H2'});
Cptrans = 2.5*R; 
dataCptrans = array2table(Cptrans);

% в расчетной работе приводятся значения и формулы для положительной
% приведенной энергии гиббса.

% расчет вращательной составляющей.
lnQrot_linear = log(T)+log(I([1,5,6],1))-log(sigma([1,3,4],1))+104.5258; 
lnQrot_nonlinear = 1.5*log(T)+0.5*log(I(2,1)*I(3,1)*I(4,1)*10^141)-log(2)-4.965; 
datalnQrotlin = array2table(lnQrot_linear,'VariableNames',{'298K' '500K' '1000K'}, 'RowNames',{'CO','CO2','H2'});
datalnQrotnonlin = array2table(lnQrot_nonlinear,'VariableNames',{'298K' '500K' '1000K'}, 'RowNames',{'H2O'});
Urot_linear = R*T; Urot_nonlinear = 3*R*T/2; 
Urot = [Urot_linear; Urot_nonlinear]; 
Hrot_linear = Urot_linear; 
Hrot_nonlinear = Urot_nonlinear; 
Hrot = [Hrot_linear;Hrot_nonlinear]; 
HH0rot_linear = R.*T;
HH0rot_nonlinear = 3*R.*T./2;
HH0rot = [HH0rot_linear; HH0rot_nonlinear];
Gibbsrot_linear = R*(log(T)+log(I([1,5,6],1))-log(sigma([1,3,4],1))+104.5258); 
Gibbsrot_nonlinear = -12.4716*log(T)+4.1572*log(I(2,1)*I(3,1)*I(4,1))-log(2)+1320.3622; 
Gibbsrot = [Gibbsrot_linear; Gibbsrot_nonlinear]; 
Srot_linear = R.*(log(T)+log(I([1,5,6],1))-log(sigma([1,3,4],1))+105.5258); 
Srot_nonlinear = R.*(1.5.*log(T)+0.5.*log(I(2,1)*I(3,1)*I(4,1))-log(sigma(2,1))+157.3621); 
Srot = [Srot_linear; Srot_nonlinear]; 
dataSrotlin = array2table(Srot_linear,'VariableNames',{'298K' '500K' '1000K'}, 'RowNames',{'CO','CO2','H2'});
dataSrotnonlin = array2table(Srot_nonlinear,'VariableNames',{'298K' '500K' '1000K'}, 'RowNames',{'H2O'});
Cprot_linear = R; 
Cprot_nonlinear = 3*R/2; 
Cprot = [Cprot_linear; Cprot_nonlinear];
dataCprot = array2table(Cprot,'VariableNames', {' '}, 'RowNames', {'Лин.', 'Нелинейн.'});


% расчет колебательной составляющей.
lnQosc_i = -log(1 - exp(-theta./T)); 
lnQosc_CO = 1*lnQosc_i(1,:); 
lnQosc_H2O = lnQosc_i(2,:)+lnQosc_i(3,:)+lnQosc_i(4,:); 
lnQosc_CO2 = lnQosc_i(5,:)+lnQosc_i(6,:)+2.*lnQosc_i(7,:); 
lnQosc_H2 = 1*lnQosc_i(8,:); 
lnQosc = [lnQosc_CO; lnQosc_H2O; lnQosc_CO2; lnQosc_H2]; 
datalnQosc = array2table(lnQosc,'VariableNames',{'298K' '500K' '1000K'}, 'RowNames',{'CO', 'H2O', 'CO2','H2'});
Uosc = R.*theta./(exp(theta./T)-1); 
Hosc = Uosc; 
HH0osc = (R.*theta)./(exp(theta./T)-1);
Cposc = (R.*(theta./T).^2.*exp(theta./T))./((exp(theta./T)-1).^2);
dataCposc = array2table(Cposc,'VariableNames',{'298K' '500K' '1000K'}, 'RowNames',{'CO', 'H2O_1', 'H2O_2', 'H2O_3', 'CO2_1', 'CO2_2', 'CO2_3','H2'});
Gibbsosc = R.*log(1-exp(-theta./T)); 
Sosc = R.*(theta./T.*1./(exp(theta./T)-1)-log(1-exp(-theta./T))); 
dataSosc = array2table(Sosc,'VariableNames',{'298K' '500K' '1000K'}, 'RowNames',{'CO', 'H2O_1', 'H2O_2', 'H2O_3', 'CO2_1', 'CO2_2', 'CO2_3','H2'});


% расчет электронной составляющей.
lnQelectron = log(1).*T.*M; 
datalnQelectron = array2table(lnQelectron,'VariableNames',{'298K' '500K' '1000K'}, 'RowNames',{'CO', 'H2O', 'CO2','H2'});
Uelectron = 0; 
Helectron = 0; 
HH0electron = 0;
Cpelectron = 0; 
Gibbselectron = -R.*log(G); 
Selectron = R.*log(1);

format shortG
% суммирование стат. сумм по каждому из веществ.
lnQ_CO = lnQtrans(1,:)+lnQrot_linear(1,:)+lnQosc(1,:)+lnQelectron(1,:); 
lnQ_H2O = lnQtrans(2,:)+lnQrot_nonlinear+lnQosc(2,:)+lnQelectron(2,:);  
lnQ_CO2 = lnQtrans(3,:)+lnQrot_linear(2,:)+lnQosc(3,:)+lnQelectron(3,:); 
lnQ_H2 = lnQtrans(4,:)+lnQrot_linear(3,:)+lnQosc(4,:)+lnQelectron(4,:);
lnQ = [lnQ_CO; lnQ_H2O; lnQ_CO2; lnQ_H2];
datalnQ = array2table(lnQ,'VariableNames',{'298K' '500K' '1000K'}, 'RowNames',{'CO', 'H2O', 'CO2','H2'});
U_CO = Utrans(1,1)+Urot(1,1)+Uosc(1,1)+Uelectron(1,1);
U_H2O = Utrans(1,1)+Urot(2,1)+Uosc(1,1)+Uelectron(1,1);
Q_CO = exp(lnQ_CO); 
Q_H2O = exp(lnQ_H2O); 
Q_CO2 = exp(lnQ_CO2); 
Q_H2 = exp(lnQ_H2);
Q = [Q_CO; Q_H2O; Q_CO2; Q_H2];

format short
% расчет единой составляющей для атомов, с более чем одной хар-ой частотой.
Uosc_H2O = Uosc(2,:)+Uosc(3,:)+Uosc(4,:);
Uosc_CO2 = Uosc(5,:)+Uosc(6,:)+2.*Uosc(7,:);
Hosc_H2O = Hosc(2,:)+Hosc(3,:)+Hosc(4,:);
Hosc_CO2 = Hosc(5,:)+Hosc(6,:)+2.*Hosc(7,:);
HH0osc_H2O = HH0osc(2,:)+HH0osc(3,:)+HH0osc(4,:);
HH0osc_CO2 = HH0osc(5,:)+HH0osc(6,:)+2.*HH0osc(7,:);
Sosc_H2O = Sosc(2,:)+Sosc(3,:)+Sosc(4,:); 
Sosc_CO2 = Sosc(5,:)+Sosc(6,:)+2.*Sosc(7,:); 
Cposc_H2O = Cposc(2,:)+Cposc(3,:)+Cposc(4,:);
Cposc_CO2 = Cposc(5,:)+Cposc(6,:)+2.*Cposc(7,:);
Gibbsosc_H2O = Gibbsosc(2,:)+Gibbsosc(3,:)+2.*Gibbsosc(4,:);
Gibbsosc_CO2 = Gibbsosc(5,:)+Gibbsosc(6,:)+2.*Gibbsosc(7,:);

% Вывод результатов: сверху вниз (CO, H2O, CO2, H2), слева направо (298, 500, 1000, 1500).
U_CO = Utrans+Urot(1,:)+Uosc(1,:)+Uelectron;
U_H2O = Utrans+Urot(2,:)+Uosc_H2O+Uelectron;
U_CO2 = Utrans+Urot(1,:)+Uosc_CO2+Uelectron;
U_H2 = Utrans+Urot(1,:)+Uosc(8,:)+Uelectron;
U = [U_CO;U_H2O;U_CO2;U_H2]; %J/mol

H_CO = Htrans+Hrot(1,:)+Hosc(1,:)+Helectron;
H_H2O = Htrans+Hrot(2,:)+Hosc_H2O+Helectron;
H_CO2 = Htrans+Hrot(1,:)+Hosc_CO2+Helectron;
H_H2 = Htrans+Hrot(1,:)+Hosc(8,:)+Helectron;
H = [H_CO;H_H2O;H_CO2;H_H2]; %J/mol

HH0_CO = HH0trans+HH0rot(1,:)+HH0osc(1,:)+HH0electron;
HH0_H2O = HH0trans+HH0rot(2,:)+HH0osc_H2O+HH0electron;
HH0_CO2 = HH0trans+HH0rot(1,:)+HH0osc_CO2+HH0electron;
HH0_H2 = HH0trans+HH0rot(1,:)+HH0osc(8,:)+HH0electron;
HH0 = [HH0_CO;HH0_H2O;HH0_CO2;HH0_H2]; %J/mol
dataHH0 = array2table(HH0,'VariableNames',{'298K' '500K' '1000K'}, 'RowNames',{'CO', 'H2O', 'CO2','H2'});

S_CO = Strans(1,:)+Srot(1,:)+Sosc(1,:)+Selectron(1,1);
S_H2O = Strans(2,:)+Srot(4,:)+Sosc_H2O+Selectron(1,1);
S_CO2 = Strans(3,:)+Srot(2,:)+Sosc_CO2+Selectron(1,1);
S_H2 = Strans(4,:)+Srot(3,:)+Sosc(8,:)+Selectron(1,1);
S = [S_CO;S_H2O;S_CO2;S_H2]; %J/mol
dataS = array2table(S,'VariableNames',{'298K' '500K' '1000K'}, 'RowNames',{'CO', 'H2O', 'CO2','H2'});


Cp_CO = Cptrans+Cprot(1,:)+Cposc(1,:)+Cpelectron;
Cp_H2O = Cptrans+Cprot(2,:)+Cposc_H2O+Cpelectron;
Cp_CO2 = Cptrans+Cprot(1,:)+Cposc_CO2+Cpelectron;
Cp_H2 = Cptrans+Cprot(1,:)+Cposc(8,:)+Cpelectron;
Cp = [Cp_CO;Cp_H2O;Cp_CO2;Cp_H2]; %J/mol
dataCp = array2table(Cp,'VariableNames',{'298K' '500K' '1000K'}, 'RowNames',{'CO', 'H2O', 'CO2','H2'});
dCp = (Cp_H2+Cp_CO2)-(Cp_CO+Cp_H2O); %J/mol


Gibbs_CO = Gibbstrans(1,:)+Gibbsrot(1,:)+Gibbsosc(1,:)+Gibbselectron(1,1);
Gibbs_H2O = R.*(lnQ(2,:)-1);
Gibbs_CO2 = Gibbstrans(3,:)+Gibbsrot(2,:)+Gibbsosc_CO2+Gibbselectron(1,1);
Gibbs_H2 = Gibbstrans(4,:)+Gibbsrot(3,:)+Gibbsosc(8,:)+Gibbselectron(1,1);
Gibbs = [Gibbs_CO;Gibbs_H2O;Gibbs_CO2;Gibbs_H2]; %J/mol
dataGibbs = array2table(Gibbs,'VariableNames',{'298K' '500K' '1000K'}, 'RowNames',{'CO', 'H2O', 'CO2','H2'});
% пересчитать для воды(155.31), для СО2(182.25), 

format short
% Расчет константы равновесия используя вычисленные значения стат. сумм.
% 1) Расчет при помощи табличных данных теплот образования.
Hform_CO = -110.53;
Hform_H2O = -241.81;
Hform_CO2 = -393.51;
Hform_H2 = 0;
Hform = ([Hform_CO; Hform_H2O; Hform_CO2; Hform_H2])*10^3;
dataHform = array2table(Hform,'RowNames',{'CO', 'H2O', 'CO2','H2'});
dHform = (Hform(4,1)+Hform(3,1))-(Hform(2,1)+Hform(1,1));

lnKa1_298 = log((Q(3,1)./exp(1)).*(Q(4,1)./exp(1))./((Q(1,1)./exp(1)).*(Q(2,1)./exp(1))))-dHform./R*298;

% 2) Расчет при помощи уравнения температурной зависимости.
dHT_298calc = (H(4,1)+H(3,1))-(H(2,1)+H(1,1));

% расчет коэффициентoв уравнения Cp = a+bT+c'/T^2
for i = [1:1:4]
    C1 = Cp(i,1);
    C2 = Cp(i,2);
    C3 = Cp(i,3);
    T1 = 298;
    T2 = 500;
    T3 = 1000;
    c(i,1) = (-((C3-C2)*(T2-T1)-(C2-C1)*(T3-T2))*(T1*T2*T3)^2)/((T2-T1)*(T3-T2)*((T3+T2)*T1^2-(T2+T1)*T3^2));
end
datac = array2table(c, 'RowNames',{'CO', 'H2O', 'CO2','H2'}, 'VariableNames',{'c`'} );

for i = [1:1:4]
        C1 = Cp(i,1);
        C2 = Cp(i,2);
        C3 = Cp(i,3);

        T1 = 298;
        T2 = 500;
        T3 = 1000;
    b(i,1) = ((C3-C2)/(T3-T2))+(c(i,1)*(T3+T2)/(T2^2*T3^2));
end
datab = array2table(b, 'RowNames',{'CO', 'H2O', 'CO2','H2'});

for i = [1:1:4]
       C1 = Cp(i,1);
        C2 = Cp(i,2);
        C3 = Cp(i,3);
        T1 = 298;
        T2 = 500;
        T3 = 1000;
    a(i,1) = C3-b(i,1)*T3-c(i,1)/T3^2;
end
dataa = array2table(a, 'RowNames',{'CO', 'H2O', 'CO2','H2'});
deltaa = ((a(4,1)+a(3,1))-(a(2,1)+a(1,1)));
deltab = ((b(4,1)+b(3,1))-(b(2,1)+b(1,1)));
deltac = ((c(4,1)+c(3,1))-(c(2,1)+c(1,1)));
deltai = [deltaa,deltab,deltac];
datadeltai = array2table(deltai, 'VariableNames',{'∆a', '∆b', '∆c`'});

%уравнение Cp = f(T):
Cp_eqCO = a(1,1)+b(1,1).*T+c(1,1)./T.^2;
Cp_eqH2O = a(2,1)+b(2,1).*T+c(2,1)./T.^2;
Cp_eqCO2 = a(3,1)+b(3,1).*T+c(3,1)./T.^2;
Cp_eqH2 = a(4,1)+b(4,1).*T+c(4,1)./T.^2;
Cp_equation = [Cp_eqCO;Cp_eqH2O;Cp_eqCO2;Cp_eqH2];
% составление зависимости dCp_equation = f(T)
% dCp_equation


%lnKa2_298 = log((Q(3,1)./exp(1)).*(Q(4,1)./exp(1))./((Q(1,1)./exp(1)).*(Q(2,1)./exp(1))))-dH0_calc./R*298

% 3) Расчет при помощи приращения энтальпии и табличного значения dH0298.
for i = [1,2,3]
dHH0(1,i) = (HH0(4,i)+HH0(3,i))-(HH0(2,i)+HH0(1,i));
end
datadHH0 = array2table(dHH0,'VariableNames',{'298K' '500K' '1000K'});

dH0 = dHform-dHH0(1,1);
dH500 = dHH0(1,2)+dH0;
dH1000 = dHH0(1,3)+dH0;



lnKa3_298 = log((Q(3,1)./exp(1)).*(Q(4,1)./exp(1))./((Q(1,1)./exp(1)).*(Q(2,1)./exp(1))))-dH0./R*298;

Ka1_298 = exp(lnKa1_298);
%Ka2_298 = exp(lnKa2_298)
Ka3_298 = exp(lnKa3_298);

% табличное значение приращения энтальпии d(HT-H0).
deltaHTH0_table = 6.20;
% kJ/mol

% расчет dH0 = dHT-d(HT-H0), dHT - расчетные значения энтальпии для 298К,
% d(HT-H0) - табличное значение из КСФХВ.
dH0_table = ((H(4,1)+H(3,1))-(H(2,1)+H(1,1)))-deltaHTH0_table*10^3;
% J/mol

dGibbs = (Gibbs(4,1)+Gibbs(3,1))-(Gibbs(2,1)+Gibbs(1,1)); %расчет изменения приведенной энергии Гиббса.

fprintf(['Таблица молекулярных постоянных:\n Молекулярная масса M, г/моль:\n CO = %.4f;\n H2O = %.4f;\n CO2 = %.4f;\n H2 = %.4f;\n Приведенная масса эквивалента m, кг:\n CO = %.4g;\n H2O = %.4g;\n CO2 = %.4g;\n H2 = %.4g;\n Межатомные растояния r, м:\n CO = %.4g;\n H2O = %.4g;\n CO2 = %.4g;\n H2 = %.4g;\n Угол в молекуле φ, радиан:\n H2O = %.4g;\n Температура T, К:\n T1 = %.f;\n T2 = %.f;\n T3 = %.f;\n Универсальная газовая постоянная, Дж/(моль*К): R = %.3f;\n Константа Больтцмана, Дж/К: k = %.4g;\n Давление, Па: P = %.f;\n Постоянная Планка, Дж*с: h = %.4g;\n'],M_CO, M_H2O, M_CO2, M_H2, m_CO, m_H2O,m_CO2, m_H2, r_CO, r_H2O, r_CO2, r_H2, angle_H2O, T,R,k,P,h) 
fprintf(['\nСоставляющие статистической суммы:\n lnQ(пост.):\n'])
disp(datalnQtrans)
writetable(datalnQtrans, 'datalnQtrans.xls')
fprintf([' lnQ(вращ. лин.):\n'])
disp(datalnQrotlin)
writetable(datalnQrotlin, 'datalnQrotlin.xls')
fprintf([' lnQ(вращ. нелин.):\n'])
disp(datalnQrotnonlin)
writetable(datalnQrotnonlin, 'datalnQrotnonlin.xls')
fprintf([' lnQ(колеб.):\n'])
disp(datalnQosc)
writetable(datalnQosc, 'datalnQosc.xls')
fprintf([' lnQ(эл.):\n'])
disp(datalnQelectron)
writetable(datalnQelectron, 'datalnQelectron.xls')

fprintf(['\nПолная статистическая сумма по состояниям:\n'])
disp(datalnQ)
writetable(datalnQ, 'datalnQ.xls')

fprintf(['\nСоставляющие энтропии, Дж/К:\n'])
fprintf([' S(пост.):\n'])
disp(dataStrans)
writetable(dataStrans, 'dataStrans.xls')
fprintf([' S(вращ. лин.):\n'])
disp(dataSrotlin)
writetable(dataSrotlin, 'dataSrotlin.xls')
fprintf([' S(вращ. нелин.):\n'])
disp(dataSrotnonlin)
writetable(dataSrotnonlin, 'dataSrotnonlin.xls')
fprintf([' S(колеб.):\n'])
disp(dataSosc)
writetable(dataSosc, 'dataSosc.xls')
fprintf([' S(эл.):\n'])
disp(Selectron)
writematrix(Selectron, 'Selectron.xls')

fprintf(['\nАбсолютная Энтропия, Дж/К\n'])
disp(dataS)
writetable(dataS, 'dataS.xls')

fprintf(['\nСоставляющие изобарной теплоемкости, Дж/(моль*К):\n'])
fprintf([' Cp(пост.):\n'])
disp(dataCptrans)
writetable(dataCptrans, 'dataCptrans.xls')
fprintf([' Cp(вращ.):\n'])
disp(dataCprot)
writetable(dataCprot, 'dataCprot.xls')
fprintf([' Cp(колеб.):\n'])
disp(dataCposc)
writetable(dataCposc, 'dataCposc.xls')
fprintf([' Cp(эл.):\n'])
disp(Cpelectron)
writematrix(Cpelectron, 'Cpelectron.xls')

fprintf(['\nИзобарная теплоемкость веществ, участвующих в реакции, Дж/(моль*К)\n'])
disp(dataCp)
writetable(dataCp, 'dataCp.xls')

fprintf(['\nВычисление коэффициентов уравнения теплоемкости Cp = f(T) и их изменения, Дж/(моль*К)\n'])
disp(dataa)
writetable(dataa, 'dataa.xls')
disp(datab)
writetable(datab, 'datab.xls')
disp(datac)
writetable(datac, 'datac.xls')
disp(datadeltai)
writetable(datadeltai, 'datadeltai.xls')
fprintf(['\nЗначения приращений энтальпии (Ht-H0), Дж/моль:\n'])
disp(dataHH0)
writetable(dataHH0, 'dataHH0.xls')
fprintf(['\nТеплоты образования веществ (табличные) ∆Hf(298), Дж/моль:\n'])
disp(dataHform)
writetable(dataHform, 'dataHform.xls')
fprintf(['\nИзменение теплот образования ∆Hf(298) = %.4f Дж/моль\n '], dHform)
fprintf(['\nИзменение приращения энтальпии ∆(Ht-H0), Дж/моль\n'])
disp(datadHH0)
writetable(datadHH0, 'datadHH0.xls')
fprintf(['\n ∆H0 = ∆Hf-∆(H298-H0) = %.4f Дж/моль\n'], dH0)
fprintf(['\nВычисленные энтальпии реакций при выбранных температурах, Дж/моль:\n ∆H(500) = %.4f;\n ∆H(1000) = %.4f;\n'], dH500, dH1000)
fprintf(['\nЗначения приращений энергии Гиббса (G0-H0)/T, Дж/моль:\n'])
disp(dataGibbs)
writetable(dataGibbs, 'dataGibbs.xls')
fprintf(['\n Изменение приращений энергии Гиббса (298K), Дж/моль: ∆(G0-H0)/T = %.4f\n'],dGibbs)
lnKa298 = -1/R*(-dGibbs+dH0/298);
A = dHform-deltaa*298-deltab/2*298^2+deltac/298;
B = lnKa298+A/(R*298)-deltaa*log(298)/R-deltab*298/(2*R)-deltac/(2*R*298^2);
fprintf(['\nРасчет термодинамических констант равновесия:\n lnKa(298) =  %.4f;\n A = %.4f;\n B = %.4f\n'], lnKa298, A, B)
fprintf(['\nРасчет термодинамических констант равновесия при других температурах:\n 1 - По уравнению lnKa1 = -1/R*(∆(G-H0)/T+∆H0/T)\n 2 - По уравнению lnKa2 = B-A/(RT)+∆a*ln(T)/R+∆bT/(2R)+∆c/(2R*T^2)\n 3 - По средним значениям теплоемкости lnKa3 = -∆G/(RT)\n'])
temp = [298:100:1000];
lnKa1T = -1/R.*(-dGibbs+dH0./temp);
lnKa2T = B-A./(R.*temp)+deltaa.*log(temp)./R+deltab./(2.*R).*temp+deltac./(2.*R.*temp.^2);
Cp_aver = [30.92; 37.06; 47.15; 29.57];
dataCp_aver = array2table(Cp_aver,'VariableNames',{'Cp(средняя),Дж/(моль*К)'},'RowNames',{'CO', 'H2O', 'CO2','H2'});
disp(dataCp_aver)
writetable(dataCp_aver, 'dataCp_aver.xls')
dCp_aver = (Cp_aver(4,1)+Cp_aver(3,1))-(Cp_aver(2,1)+Cp_aver(1,1));
fprintf(['\nИзменение средней теплоемкости ∆Cp(средн.) = %.4f Дж/(моль*К);\n'], dCp_aver)
dH3 = dHform+dCp_aver.*(temp-298);
dS3 = ((S(4,1)+S(3,1))-(S(2,1)+S(1,1)))+dCp_aver.*log(temp./298);
dG3 = dH3-temp.*dS3;
lnKa3T = -dG3./(R.*temp);
Constants = transpose([[300:100:1000]; lnKa1T;lnKa2T;lnKa3T]);
datalnKa = array2table(Constants, 'VariableNames',{'T, K', 'lnKa1', 'lnKa2','lnKa3'});
disp(datalnKa)
writetable(datalnKa, 'datalnKa.xls')

figure
plot(temp, lnKa1T, '-.+'), grid on, axis padded;
title('Зависимость константы равновесия от температуры')
hold on
plot(temp, lnKa2T, '--o'), grid on, axis padded
hold on
plot(temp, lnKa3T, '-*'), grid on, axis padded
xlabel('T, K')
ylabel('lnK_{a}')
legend('lnK_{a1}','lnK_{a2}','lnK_{a3}', 'Location','best')
hold off
syms ksi Ka1
eqn1 =  Ka1 == ksi^2/(1-ksi)^2;
Ka1 = exp(lnKa1T);
S =  solve(eqn1,ksi,'ReturnConditions',true);
Pstandart = 1;
nu1 =1; nu2 = 1; nu3 = 1; nu4 = 1;
dnu = 0;
n01 = 1; n02 = 1; n03 = 0; n04 = 0;
n0i = [n01; n02; n03; n04];
fprintf(['\nРасчет химических равновесий\n\n  Исходный состав n0i, моль:\n CO = 1;\n H2O = 1;\n CO2 = 0;\n H2 = 0;\n\n  Равновесный состав ni, моль:\n CO = n0(CO)-ν(CO)* ξ;\n H2O = n0(H2O)-ν(H2O)* ξ;\n CO2 = n0(CO2)+ν(CO2)* ξ;\n H2 = n0(H2)+ν(H2)* ξ;\n ∆v, моль = %.4f\n '], dnu)
fprintf(['\n Расчет исходя из уравнений:\n\n '])
disp(eqn1)
disp(['ξ ≤ 1'])
disp('ξ = ')
disp(S.ksi)
ksi1 = Ka1.^(1/2)./(Ka1.^(1/2) - 1);
ksi2 = Ka1.^(1/2)./(Ka1.^(1/2) + 1);
ksi_table = [ksi2(1,1) ksi2(1,3) ksi2(1,8)];
fprintf(['\n При 1 атм. и значениях Ka1(T), ξ принимает значения:\n'])
dataksi = array2table(ksi_table,'VariableNames',{'298K' '500K' '1000K'}, 'RowNames', {'ξ'});
disp(dataksi)
writetable(dataksi, 'dataksi.xls')
fprintf(['\nСостав равновесной смеси в молях\n T = 298K:\n CO = %d-%.4f;\n H2O = %d-%.4f;\n CO2 = %d+%.4f;\n H2 = %d+%.4f;\n\n T = 500K:\n CO = %d-%.4f;\n H2O = %d-%.4f;\n CO2 = %d+%.4f;\n H2 = %d+%.4f;\n\n T = 1000K:\n CO = %d-%.4f;\n H2O = %d-%.4f;\n CO2 = %d+%.4f;\n H2 = %d+%.4f;\n'], n01,ksi2(1,1), n02,ksi2(1,1), n03,ksi2(1,1), n04,ksi2(1,1), n01,ksi2(1,1), n02,ksi2(1,3), n03,ksi2(1,3), n04,ksi2(1,3),n01,ksi2(1,8), n02,ksi2(1,8), n03,ksi2(1,8), n04,ksi2(1,8))
fprintf(['\nСостав равновесной смеси в мольных долях\nT = 298K:\n'])
sumni298 = (n01-ksi2(1,1))+(n02-ksi2(1,1))+(n03+ksi2(1,1))+(n04+ksi2(1,1));
Xi298 = [(n01-ksi2(1,1))/sumni298; (n02-ksi2(1,1))/sumni298; (n03+ksi2(1,1))/sumni298; (n04+ksi2(1,1))/sumni298];
dataXi298 = array2table(Xi298, 'RowNames', {'CO','H2O','CO2', 'H2'});
disp(dataXi298)
writetable(dataXi298, 'dataXi298.xls')
fprintf(['\nT = 500K:\n'])
sumni500 = (n01-ksi2(1,3))+(n02-ksi2(1,3))+(n03+ksi2(1,3))+(n04+ksi2(1,3));
Xi500 = [(n01-ksi2(1,3))/sumni500; (n02-ksi2(1,3))/sumni500; (n03+ksi2(1,3))/sumni500; (n04+ksi2(1,3))/sumni500];
dataXi500 = array2table(Xi500, 'RowNames', {'CO','H2O','CO2', 'H2'});
disp(dataXi500)
writetable(dataXi500, 'dataXi500.xls')
fprintf(['\nT = 1000K:\n'])
sumni1000 = (n01-ksi2(1,8))+(n02-ksi2(1,8))+(n03+ksi2(1,8))+(n04+ksi2(1,8));
Xi1000 = [(n01-ksi2(1,8))/sumni1000; (n02-ksi2(1,8))/sumni1000; (n03+ksi2(1,8))/sumni1000; (n04+ksi2(1,8))/sumni1000];
dataXi1000 = array2table(Xi1000, 'RowNames', {'CO','H2O','CO2', 'H2'});
disp(dataXi1000)
writetable(dataXi1000, 'dataXi1000.xls')
fprintf(['\nПарциальное давление компонентов смеси при T = 298K:\n'])
pi = Xi298.*Pstandart;
datapi = array2table(pi,'RowNames', {'CO','H2O','CO2', 'H2'});
disp(datapi)
writetable(datapi, 'datapi.xls')
alpha = ksi2(1,1)/n01;
fprintf(['\nСтепень превращения  исходного вещества при T = 298K: α = %.4f\n'], alpha)
beta = (n03+ksi2(1,1))/sumni298;
fprintf(['\nВыход продукта реакции при T = 298K: β = %.4f\n'],beta)
fprintf(['\n\n Формирование кристаллического зародыша.\n Термодинамические данные, полученные при помощи программы ИвтанТермо:\n'])
IvtanThermo = [  -580.2434e3         37.1976       -585.7860e3        -18.4758e3        101.0274;
                 -578.6895e3          9.1418       -584.2451e3        -13.8888e3         75.5676;
                 -577.3849e3          3.5478       -583.6190e3        -12.4684e3         60.3178;
                 -576.1566e3         -1.0042       -583.5004e3        -12.2399e3         50.1579;
                 -574.9150e3         -4.5396       -583.7858e3        -12.6727e3         42.8999;
                 -573.6105e3         -7.1596       -584.3778e3        -13.4595e3         37.4522;
                 -572.2175e3         -8.9939       -585.1916e3        -14.4156e3         33.2100;
                 -570.7253e3        -10.1844       -586.1553e3        -15.4299e3         29.8111];
dataIvtanThermo = array2table(IvtanThermo, 'VariableNames', {'∆G(T) Дж/моль','∆Cp(T) кДж/(моль*K)','∆H(T) Дж/моль','∆S(T) Дж/(моль*K)','log10(Kp)'},'RowNames', {'300','400','500','600','700','800','900','1000'});
disp(dataIvtanThermo)
writetable(dataIvtanThermo, 'dataIvtanThermo.xls')
figure
plot([300:100:1000],IvtanThermo(:,1)), grid on, axis padded;
xlabel('T, K')
ylabel('∆G(T), Дж/моль')
legend('∆G = f(T)', 'Location','best')
sigma300 = 10.07; %J
fprintf(['\nТабличное значение поверхностного натяжения J = %.2f Н/м\n'],sigma300)

sigma400 = sigma300-0.02.*sigma300;
sigma500 = sigma400-0.02.*sigma400;
sigma600 = sigma500-0.02.*sigma500;
sigma700 = sigma600-0.02.*sigma600;
sigma800 = sigma700-0.02.*sigma700;
sigma900 = sigma800-0.02.*sigma800;
sigma1000 = sigma900-0.02.*sigma900;

sigmaT = [sigma400;sigma500;sigma600;sigma700;sigma800;sigma900;sigma1000]; %J
datasigmaT = array2table(sigmaT, 'RowNames', {'400', '500', '600', '700','800','900','1000'});
fprintf(['\nЗначение показателя поверхностного натяжения при заданных температурах:\n'])
disp(datasigmaT)
writetable(datasigmaT, 'datasigmaT.xls')
figure
plot([400:100:1000],sigmaT), grid on, axis padded
xlabel('T, K')
ylabel('σ, Дж')
legend('σ = f(T)', 'Location','best')
fprintf(['\nВлияние размера кристаллического зародыша на поверхностную, объемную и суммарную составляющие энергии Гиббса процесса зародышеобразования. Определение критического размера зародыша r(кр) и максимальной энергии Гиббса ∆Gmax\n Примем, что зародыш имеет форму сферы.\n'])
r = [1:1:15]*10^-9; %m
t = [400:100:1000]; %K
gamma = 4; %sphere
dGs = 3.1415926*gamma.*transpose(r).^2.*(transpose(sigmaT)./8);
rho_Cr2O3 = 5.22; %g/cm3
M_Cr2O3 = 151.9904;%g/mol
Vcr = 4/3*3.1415926.*r.^3; %m3
Vm = rho_Cr2O3/M_Cr2O3; %cm3 далее переведено внутри формул в м3
Tmelt_Cr2O3 = 2708; %K
dHmelt_Cr2O3 = 125000;%J/mol
dGv = -(transpose(Vcr).*rho_Cr2O3*10^6.*(Tmelt_Cr2O3-t./Tmelt_Cr2O3));
dGsum = dGv+dGs;

fprintf(['\nРезультаты расчета ∆Gs, ∆Gv, ∆Gsum в зависимости от размеров зародыша при заданной температуре представленны в следующих таблицах:\n'])

fprintf(['\nT = 400K;\n'])
CrystGibbs400 = [transpose(r) dGs(:,1) dGv(:,1) dGsum(:,1)];
dataCrystGibbs400 = array2table(CrystGibbs400,'VariableNames', {'r, m','∆Gs, Дж/моль','∆Gv, Дж/моль','∆Gsum, Дж/моль'});
disp(dataCrystGibbs400)
writetable(dataCrystGibbs400, 'dataCrystGibbs400.xls')
fprintf(['\nT = 500K;\n'])
CrystGibbs500 = [transpose(r) dGs(:,2) dGv(:,2) dGsum(:,2)];
dataCrystGibbs500 = array2table(CrystGibbs500,'VariableNames', {'r, m','∆Gs, Дж/моль','∆Gv, Дж/моль','∆Gsum, Дж/моль'});
disp(dataCrystGibbs500)
writetable(dataCrystGibbs500, 'dataCrystGibbs500.xls')
fprintf(['\nT = 600K;\n'])
CrystGibbs600 = [transpose(r) dGs(:,3) dGv(:,3) dGsum(:,3)];
dataCrystGibbs600 = array2table(CrystGibbs600,'VariableNames', {'r, m','∆Gs, Дж/моль','∆Gv, Дж/моль','∆Gsum, Дж/моль'});
disp(dataCrystGibbs600)
writetable(dataCrystGibbs600, 'dataCrystGibbs600.xls')
fprintf(['\nT = 700K;\n'])
CrystGibbs700 = [transpose(r) dGs(:,4) dGv(:,4) dGsum(:,4)];
dataCrystGibbs700 = array2table(CrystGibbs700,'VariableNames', {'r, m','∆Gs, Дж/моль','∆Gv, Дж/моль','∆Gsum, Дж/моль'});
disp(dataCrystGibbs700)
writetable(dataCrystGibbs700, 'dataCrystGibbs700.xls')
fprintf(['\nT = 800K;\n'])
CrystGibbs800 = [transpose(r) dGs(:,5) dGv(:,5) dGsum(:,5)];
dataCrystGibbs800 = array2table(CrystGibbs800,'VariableNames', {'r, m','∆Gs, Дж/моль','∆Gv, Дж/моль','∆Gsum, Дж/моль'});
disp(dataCrystGibbs800)
writetable(dataCrystGibbs800, 'dataCrystGibbs800.xls')
fprintf(['\nT = 900K;\n'])
CrystGibbs900 = [transpose(r) dGs(:,6) dGv(:,6) dGsum(:,6)];
dataCrystGibbs900 = array2table(CrystGibbs900,'VariableNames', {'r, m','∆Gs, Дж/моль','∆Gv, Дж/моль','∆Gsum, Дж/моль'});
disp(dataCrystGibbs900)
writetable(dataCrystGibbs900, 'dataCrystGibbs900.xls')
fprintf(['\nT = 1000K;\n'])
CrystGibbs1000 = [transpose(r) dGs(:,7) dGv(:,7) dGsum(:,7)];
dataCrystGibbs1000 = array2table(CrystGibbs1000,'VariableNames', {'r, m','∆Gs, Дж/моль','∆Gv, Дж/моль','∆Gsum, Дж/моль'});
disp(dataCrystGibbs1000)
writetable(dataCrystGibbs1000, 'dataCrystGibbs1000.xls')


figure
subplot(2,4,1)
plot(r,dGs(:,1)), grid on, axis padded
xlabel('r, m')
ylabel('∆G_i, T = 400, Дж/моль')
hold on
plot(r,dGv(:,1)), grid on, axis padded
hold on
plot(r,dGsum(:,1)), grid on, axis padded
legend('∆G_s = f(r)','∆G_v = f(r)','∆G_{сумм} = f(r)', 'Location','best')
hold off

subplot(2,4,2)
plot(r,dGs(:,2)), grid on, axis padded
xlabel('r, m')
ylabel('∆G_i, T = 500, Дж/моль')
hold on
plot(r,dGv(:,2)), grid on, axis padded
hold on
plot(r,dGsum(:,2)), grid on, axis padded
legend('∆G_s = f(r)','∆G_v = f(r)','∆G_{сумм} = f(r)', 'Location','best')
hold off

subplot(2,4,3)
plot(r,dGs(:,3)), grid on, axis padded
xlabel('r, m')
ylabel('∆G_i, T = 600, Дж/моль')
hold on
plot(r,dGv(:,3)), grid on, axis padded
hold on
plot(r,dGsum(:,3)), grid on, axis padded
legend('∆G_s = f(r)','∆G_v = f(r)','∆G_{сумм} = f(r)', 'Location','best')
hold off

subplot(2,4,4)
plot(r,dGs(:,4)), grid on, axis padded
xlabel('r, m')
ylabel('∆G_i, T = 700, Дж/моль')
hold on
plot(r,dGv(:,4)), grid on, axis padded
hold on
plot(r,dGsum(:,4)), grid on, axis padded
legend('∆G_s = f(r)','∆G_v = f(r)','∆G_{сумм} = f(r)', 'Location','best')
hold off

subplot(2,4,5)
plot(r,dGs(:,5)), grid on, axis padded
xlabel('r, m')
ylabel('∆G_i, T = 800, Дж/моль')
hold on
plot(r,dGv(:,5)), grid on, axis padded
hold on
plot(r,dGsum(:,5)), grid on, axis padded
legend('∆G_s = f(r)','∆G_v = f(r)','∆G_{сумм} = f(r)', 'Location','best')
hold off

subplot(2,4,6)
plot(r,dGs(:,6)), grid on, axis padded
xlabel('r, m')
ylabel('∆G_i, T = 900, Дж/моль')
hold on
plot(r,dGv(:,6)), grid on, axis padded
hold on
plot(r,dGsum(:,6)), grid on, axis padded
legend('∆G_s = f(r)','∆G_v = f(r)','∆G_{сумм} = f(r)', 'Location','best')
hold off

subplot(2,4,7)
plot(r,dGs(:,7)), grid on, axis padded
xlabel('r, m')
ylabel('∆G_i, T = 1000, Дж/моль')
hold on
plot(r,dGv(:,7)), grid on, axis padded
hold on
plot(r,dGsum(:,7)), grid on, axis padded
legend('∆G_s = f(r)','∆G_v = f(r)','∆G_{сумм} = f(r)', 'Location','best')
hold off


fprintf(['\nВлияние температуры на размер критического зародыша:\n'])
rcrit =  2.*sigmaT.*Vm./-IvtanThermo([2,3,4,5,6,7,8],1);
datarcrit = array2table(rcrit,"RowNames", {'400', '500', '600', '700','800','900','1000'}, 'VariableNames',{'r(кр.)'});
disp(datarcrit)

% fprintf(['\n\n'])