clc; clear all;
Isp1 = 250;
Isp2 = 300;
g = 9.81;
R_earth = 6553.630415; % km
miu_earth = 398600.4; 
Vbo = 1000*sqrt(miu_earth / R_earth); %burnout velocity 
ml = 7.598857008462339e+03;   %chosen payload mass
sig1= .07;  %structure ratio first stage     sig = mst/(mst+mpr)
sig2 = .05; %structure ratio second stage



r = 0:.01:1;
syms 'lam1' 'lam2';

lameqn = Vbo+(g*Isp1*log(sig1+(1-sig1)*lam1)+(g*Isp2*log(sig2+(1-sig2)*lam2)));

lamfnc = matlabFunction(solve(lameqn,lam2));

lamtmax = max(lamfnc(r).*r);
rmax = r(lamfnc(r).*r == lamtmax);

op = plot(r,lamfnc(r),r,lamfnc(r).*r,rmax,lamtmax,'ok', rmax, lamfnc(rmax),'ok');

legend('\lambda_2','\lambda_t');
xlabel('\lambda_1');
ylabel('\lambda_2, \lambda_t');
title('Two Stage Payload Ratio Optimization');

txt = ['\lambda_1: ' num2str(rmax) ];
text(rmax-.3,lamfnc(rmax)+.5,txt);

txt2 = ['\lambda_2: ' num2str(lamfnc(rmax)) ];
text(rmax-.3,lamfnc(rmax)+.4,txt2);

txt3 = ['\lambda_t: ' num2str(lamfnc(rmax)*rmax) ];
text(rmax-.3,lamfnc(rmax)+.3,txt3);

lamt = lamtmax;
lam1 = rmax;
lam2 = lamt/rmax;

%% mass calculations
%ms# = structure mass
%mp# = propellant mass

syms ms1 ms2 mp1 mp2
eqn1 = 0 == (ms2+mp2+ml)/(mp1+ms1+mp2+ms2+ml) - lam1;
eqn2 = 0 == ml/(mp2+ms2+ml) - lam2;
eqn4 = 0 == ms1/(ms1+mp1) - sig1;
eqn5 = 0 == ms2/(ms2+mp2) - sig2;

a = solve([eqn1 eqn2 eqn4 eqn5], [ms1 ms2 mp1 mp2]);

ms1 = double(a.ms1);
mp1 = double(a.mp1);

ms2 = double(a.ms2);
mp2 = double(a.mp2);
 
 
launchMass = mp1+ms1+mp2+ms2+ml;
 
secondStageMass = ms2+mp2+ml;

 
txt4 = ['Payload Mass: ' num2str(ml) ' kg' newline 'Launch Mass: ' num2str(round(launchMass,2)) ' kg' newline 'Second Stage Mass: ' num2str(round(secondStageMass,2)) ' kg' newline ];
text(rmax-.1,lamfnc(rmax)+.5,txt4);

txt5 = ['ms1: ' num2str(ms1) ' kg; mp1: ' num2str(mp1) ' kg' newline 'ms2: ' num2str(ms2) ' kg; mp2: ' num2str(mp2) ' kg' newline ];
text(rmax-.1,lamfnc(rmax)+.2,txt5);