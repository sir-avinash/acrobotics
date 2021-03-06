%% acrobotics_main

% clear all
close all
addpath('/home/exx/Avinash/My Toolbox');
force_save = 1;


%%%%%%%%%%%%% System Properties %%%%%%%%%%
data.g = 9.81;
data.e3 = [0 0 1]';
data.m = 1;
data.mt = 0;
data.l = 1;
data.J = 1e-3*diag([0,0,1]); %eye(3);                % pendulum inertia w.r.t body-frame
data.theta1 = deg2rad(-10);     % Pitch := Rotation about Y-axis; Start with some negative deflection
data.phi1 = deg2rad(0);        % Roll := Rotation about X-axis; Start with some positive displacement
data.theta2 = deg2rad(-10);     % Pitch := Rotation about Y-axis; Start with some negative deflection
data.phi2 = deg2rad(0);        % Roll := Rotation about X-axis; Start with some positive displacement
data.dot_theta1 = deg2rad(-3);  
data.dot_phi1 = deg2rad(0);
data.dot_theta2 = deg2rad(-3);
data.dot_phi2 = deg2rad(0);

flag.Rw2b = 1; % convert omega to dph dth dps

%%%%%%%%%% Defining Initial Conditions Euler : Defining from world-frame to base-frame %%%%%%%%%
th1 = data.theta1;
ph1 = data.phi1;
th2 = data.theta2; 
ph2 = data.phi2;

dph1 = data.dot_phi1; 
dth1 = data.dot_theta1;
dph2 = data.dot_phi2;
dth2 = data.dot_theta2;

x0_eu = [ph1;th1;ph2;th2;dph1;dth1;dph2;dth2];

%%%%%%%%%% Defining Initial Conditions S^2 : Defining from world-frame to base-frame %%%%%%%%%
q10 = Rx(data.phi1)*Ry(data.theta1)*data.e3;
q20 = Rx(data.phi2)*Ry(data.theta2)*data.e3;

dq10 = [dth1*cos(th1);dth1*sin(ph1)*sin(th1) - dph1*cos(ph1)*cos(th1);- dph1*cos(th1)*sin(ph1) - dth1*cos(ph1)*sin(th1)];
dq20 = [dth2*cos(th2);dth2*sin(ph2)*sin(th2) - dph2*cos(ph2)*cos(th2);- dph2*cos(th2)*sin(ph2) - dth2*cos(ph2)*sin(th2)];
   

if(flag.Rw2b)       
w10 = hat(q10)*dq10; w20 = hat(q20)*dq20;
else
w10 = [dph1;dth1;0]; w20 = [dph2;dth2;0];
end
%            
x0_s2 = [q10;q20;w10;w20];

%%%%%%%%%% Defining Initial Conditions SO(3) : Defining from world-frame to base-frame %%%%%%%%%
R10 = Rx(data.phi1)*Ry(data.theta1);
R20 = Rx(data.phi2)*Ry(data.theta2);

Om10 = [dph1;0;0] + Rx(ph1)*[0;dth1;0]; 
Om20 = [dph2;0;0] + Rx(ph2)*[0;dth2;0]; 
           
x0_so3 = [reshape(R10,9,1);Om10;reshape(R20,9,1);Om20];


%%%%%%%%%% ODE Setup %%%%%%%%% 

options = odeset('RelTol',1e-7,'AbsTol',1e-8); 
Tend = 1;
[T_s2,X_s2]=ode45(@acrobot_sim_s2,[0 Tend],x0_s2,options,data);
[T_so3,X_so3]=ode45(@acrobot_sim_so3,[0 Tend],x0_so3,options,data);
[T_eu,X_eu]=ode45(@acrobot_sim_euler,[0 Tend],x0_eu,options,data);


Fs=80;
[T_s2,X_s2] = even_sample(T_s2,X_s2,Fs);
[T_so3,X_so3] = even_sample(T_so3,X_so3,Fs);
[T_eu,X_eu] = even_sample(T_eu,X_eu,Fs);

%%%%%%%%%% Saving data for plotting %%%%%%%%% 
if length(T_s2) == length(T_so3) && length(T_s2) == length(T_eu)
    T = T_s2;
else
    error('Times dont match after even sampling!')
end

for i=1:length(T)
    
%% Logging spherical coordinate data    
spherical.ph1(i) = atan2(-X_s2(i,2),X_s2(i,3));
spherical.th1(i) = atan2(X_s2(i,1),sqrt(X_s2(i,2)^2+X_s2(i,3)^2)); 
spherical.ph2(i) = atan2(-X_s2(i,5),X_s2(i,6));
spherical.th2(i) = atan2(X_s2(i,4),sqrt(X_s2(i,5)^2+X_s2(i,6)^2));

spherical.q1(:,i) = X_s2(i,1:3)';
spherical.q2(:,i) = X_s2(i,4:6)';
spherical.w1(:,i) = X_s2(i,7:9)';
spherical.w2(:,i) = X_s2(i,10:12)';

if(flag.Rw2b)       

spherical.dq1(:,i) = hat(spherical.w1(:,i))*spherical.q1(:,i); 
spherical.dq2(:,i) = hat(spherical.w2(:,i))*spherical.q2(:,i);
spherical.dth1(i) = spherical.dq1(1,i)/cos(spherical.th1(i));
spherical.dph1(i) = (-(spherical.dq1(2,i) + spherical.dq1(3,i)) + spherical.dth1(i)*sin(spherical.th1(i))*(sin(spherical.ph1(i))-cos(spherical.ph1(i))))/(cos(spherical.th1(i))*(sin(spherical.ph1(i))+cos(spherical.ph1(i))));
spherical.dth2(i) = spherical.dq2(1,i)/cos(spherical.th2(i));
spherical.dph2(i) = (-(spherical.dq2(2,i) + spherical.dq2(3,i)) + spherical.dth2(i)*sin(spherical.th2(i))*(sin(spherical.ph2(i))-cos(spherical.ph2(i))))/(cos(spherical.th2(i))*(sin(spherical.ph2(i))+cos(spherical.ph2(i))));

else
dph1(i) = w1(1,i);
dth1(i) = w1(2,i);
dph2(i) = w2(1,i);
dth2(i) = w2(2,i);
end

%% Logging euler coordinate data
euler.ph1(i) = wrapToPi(X_eu(i,1));
euler.th1(i) = wrapToPi(X_eu(i,2));
euler.ph2(i) = wrapToPi(X_eu(i,3));
euler.th2(i) = wrapToPi(X_eu(i,4));
euler.dph1(i) = X_eu(i,5);
euler.dth1(i) = X_eu(i,6);
euler.dph2(i) = X_eu(i,7);
euler.dth2(i) = X_eu(i,8);
euler.q1(:,i) = [sin(euler.th1(i));-cos(euler.th1(i))*sin(euler.ph1(i));cos(euler.th1(i))*cos(euler.ph1(i))];
euler.q2(:,i) = [sin(euler.th2(i));-cos(euler.th2(i))*sin(euler.ph2(i));cos(euler.th2(i))*cos(euler.ph2(i))];
euler.dq1(:,i) = [cos(euler.th1(i))*euler.dth1(i);sin(euler.ph1(i))*sin(euler.th1(i))*euler.dth1(i)-cos(euler.ph1(i))*cos(euler.th1(i))*euler.dph1(i);-(sin(euler.ph1(i))*cos(euler.th1(i))*euler.dph1(i)+cos(euler.ph1(i))*sin(euler.th1(i))*euler.dth1(i))];
euler.dq2(:,i) = [cos(euler.th2(i))*euler.dth2(i);sin(euler.ph2(i))*sin(euler.th2(i))*euler.dth2(i)-cos(euler.ph2(i))*cos(euler.th2(i))*euler.dph2(i);-(sin(euler.ph2(i))*cos(euler.th2(i))*euler.dph2(i)+cos(euler.ph2(i))*sin(euler.th2(i))*euler.dth2(i))];
euler.w1(:,i) = hat(euler.q1(:,i))*euler.dq1(:,i);
euler.w2(:,i) = hat(euler.q2(:,i))*euler.dq2(:,i);
X_eu_updated(i,:) = [euler.q1(:,i)' euler.q2(:,i)'];

%% Logging special orthognal coordinate data

special.R1(:,:,i) = reshape(X_so3(i,1:9),3,3);
special.Om1(:,i) = X_so3(i,10:12)';
special.R2(:,:,i) = reshape(X_so3(i,13:21),3,3);
special.Om2(:,i) = X_so3(i,22:24)';

[special.ph1(i),special.th1(i),~] = Rot2Eul(special.R1(:,:,i),'radians');
[special.ph2(i),special.th2(i),~] = Rot2Eul(special.R2(:,:,i),'radians');

special.q1(:,i) = special.R1(:,:,i)*data.e3;
special.q2(:,i) = special.R2(:,:,i)*data.e3;

special.dq1(:,i) = special.R1(:,:,i)*hat(special.Om1(:,i))*data.e3;
special.dq2(:,i) = special.R2(:,:,i)*hat(special.Om2(:,i))*data.e3;

X_so3_updated(i,:) = [special.q1(:,i)' special.q2(:,i)'];


% error.ph1(i) = euler.ph1(i) - ph1(i);
% error.th1(i) = euler.th1(i) - th1(i);
% error.ph2(i) = euler.ph2(i) - ph2(i);
% error.th2(i) = euler.th2(i) - th2(i);
% error.dph1(i) = euler.dph1(i) - dph1(i);
% error.dth1(i) = euler.dth1(i) - dth1(i);
% error.dph2(i) = euler.dph2(i) - dph2(i);
% error.dth2(i) = euler.dth2(i) - dth2(i);


%%%%%% Computing Energies %%%%%%%%%%%%%%%%
m1 = data.m;
m2 = data.m;
l1 = data.l;
l2 = data.l;
g  = data.g;
e3 = data.e3;

%% Energy in S^2
spherical.KE(i) = 0.5*(m1*(l1^2)*(spherical.dq1(:,i)'*spherical.dq1(:,i)) + ...
        m2*(l1*spherical.dq1(:,i)+l2*spherical.dq2(:,i))'*(l1*spherical.dq1(:,i)+l2*spherical.dq2(:,i))); 
spherical.PE(i) = m1*g*l1*spherical.q1(:,i)'*e3 + m2*g*(l1*spherical.q1(:,i)+l2*spherical.q2(:,i))'*e3 ;
spherical.TE(i) = spherical.KE(i) + spherical.PE(i) ;

%% Energy in Euler
euler.KE(i) = 0.5*(m1*(l1^2)*(euler.dq1(:,i)'*euler.dq1(:,i)) + ...
        m2*(l1*euler.dq1(:,i)+ l2*euler.dq2(:,i))'*(l1*euler.dq1(:,i)+ l2*euler.dq2(:,i))); 
euler.PE(i) = m1*g*l1*euler.q1(:,i)'*e3 + m2*g*(l1*euler.q1(:,i)+l2*euler.q2(:,i))'*e3 ;
euler.TE(i) = euler.KE(i) + euler.PE(i) ;

% %% Alternate Energy calculation in Euler
% euler2.KE(i) = acrobot_kinetic_energy(euler.dph1(i),euler.dph2(i),euler.dth1(i),euler.dth2(i),l1,l2,m1,m2,euler.ph1(i),euler.ph2(i),euler.th1(i),euler.th2(i));
% euler2.PE(i) = acrobot_potential_energy(g,l1,l2,m1,m2,euler.ph1(i),euler.ph2(i),euler.th1(i),euler.th2(i));
% euler2.TE(i) = euler2.KE(i) + euler2.PE(i);
%% Energy in SO(3)
special.KE(i) = 0.5*(m1*(l1^2)*(special.dq1(:,i)'*special.dq1(:,i)) + ...
        m2*(l1*special.dq1(:,i)+ l2*special.dq2(:,i))'*(l1*special.dq1(:,i)+ l2*special.dq2(:,i))); 
special.PE(i) = m1*g*l1*special.q1(:,i)'*e3 + m2*g*(l1*special.q1(:,i) + l2*special.q2(:,i))'*e3 ;
special.TE(i) = special.KE(i) + special.PE(i) ;

disp(['Time step no: ',num2str(i)])
end

if exist('main_data')==1 && force_save==1
save('main_data')
end

%%%% Plots 

figure(1)
subplot(3,1,1)
plot(T,euler.KE, T,spherical.KE, T,special.KE)
title('Kinetic Energy')
subplot(3,1,2)
plot(T,euler.PE, T,spherical.PE, T,special.PE)
title('Potential Energy')
subplot(3,1,3)
plot(T,euler.TE, T,spherical.TE, T,special.TE)
title('Total Energy')
legend('Euler','S^2','SO(3)')

figure(2)
subplot(2,2,1)
plot(T,rad2deg(euler.ph1), T,rad2deg(spherical.ph1), T,rad2deg(special.ph1))
title('\phi_1')
subplot(2,2,2)
plot(T,rad2deg(euler.ph2), T,rad2deg(spherical.ph2), T,rad2deg(special.ph2))
title('\phi_2')
subplot(2,2,3)
plot(T,rad2deg(euler.th1), T,rad2deg(spherical.th1), T,rad2deg(special.th1))
title('\theta_1')
subplot(2,2,4)
plot(T,rad2deg(euler.th2), T,rad2deg(spherical.th2), T,rad2deg(special.th2))
title('\theta_2')
legend('Euler','S^2','SO(3)')

figure(4)
subplot(4,1,1)
plot(T,euler.q1(1,:), T,spherical.q1(1,:), T,special.q1(1,:))
title('q1(1)')
subplot(4,1,2)
plot(T,euler.q1(2,:), T,spherical.q1(2,:), T,special.q1(2,:))
title('q1(2)')
subplot(4,1,3)
plot(T,euler.q1(3,:), T,spherical.q1(3,:), T,special.q1(3,:))
title('q1(3)')
subplot(4,1,4)
plot(T,vecnorm(euler.q1), T,vecnorm(spherical.q1), T,vecnorm(special.q1))
title('norm(q1)')
legend('Euler','S^2','SO(3)')

figure(5)
subplot(4,1,1)
plot(T,euler.q2(1,:), T,spherical.q2(1,:), T,special.q2(1,:))
title('q2(1)')
subplot(4,1,2)
plot(T,euler.q2(2,:), T,spherical.q2(2,:), T,special.q2(2,:))
title('q2(2)')
subplot(4,1,3)
plot(T,euler.q2(3,:), T,spherical.q2(3,:), T,special.q2(3,:))
title('q2(3)')
subplot(4,1,4)
plot(T,vecnorm(euler.q2), T,vecnorm(spherical.q2), T,vecnorm(special.q2))
title('norm(q2)')
legend('Euler','S^2','SO(3)')

figure(6)
subplot(3,1,1)
plot(T,euler.w1(1,:), T,spherical.w1(1,:), T,special.Om1(1,:))
title('w1(1)')
subplot(3,1,2)
plot(T,euler.w1(2,:), T,spherical.w1(2,:), T,special.Om1(2,:))
title('w1(2)')
subplot(3,1,3)
plot(T,euler.w1(3,:), T,spherical.w1(3,:), T,special.Om1(3,:))
title('w1(3)')
legend('Euler','S^2','SO(3)')

figure(7)
subplot(3,1,1)
plot(T,euler.w2(1,:), T,spherical.w2(1,:), T,special.Om2(1,:))
title('w2(1)')
subplot(3,1,2)
plot(T,euler.w2(2,:), T,spherical.w2(2,:), T,special.Om2(2,:))
title('w2(2)')
subplot(3,1,3)
plot(T,euler.w2(3,:), T,spherical.w2(3,:), T,special.Om2(3,:))
title('w2(3)')
legend('Euler','S^2','SO(3)')

% figure(3)
% subplot(2,1,1)
% plot(T,spherical.dph1,T,euler.dph1); grid on;
% title('\phi_1 Vel')
% subplot(2,1,2)
% plot(T,spherical.dph2,T,euler.dph2) ; grid on;
% title('\phi_2 Vel')
% legend('S2','Euler')
% 
% figure(4)
% subplot(2,1,1)
% plot(T,spherical.dth1,T,euler.dth1) ; grid on;
% title('\theta_1 Vel')
% subplot(2,1,2)
% plot(T,spherical.dth2,T,euler.dth2) ; grid on;
% title('\theta_2 Vel')
% legend('S2','Euler')
% 

% figure(5)
% subplot(2,1,1)
% plot(T,error.th1,T,error.ph1,T,error.th2,T,error.ph2)
% legend('e_{\theta_1}','e_{\phi_1}','e_{\theta_2}','e_{\phi_2}')
% subplot(2,1,2)
% plot(T,error.dth1,T,error.dph1,T,error.dth2,T,error.dph2)
% legend('e_{\theta_1} dot','e_{\phi_1} dot','e_{\theta_2} dot','e_{\phi_2} dot')

% figure
% subplot(3,1,1)
% plot(T,Om1(1,:),T,euler.Om1(1,:))
% title('\Omega_{1x}')
% subplot(3,1,2)
% plot(T,Om1(2,:),T,euler.Om1(2,:))
% title('\Omega_{1y}')
% subplot(3,1,3)
% plot(T,Om1(3,:),T,euler.Om1(3,:))
% title('\Omega_{1z}')
% legend('S2','Euler')
% 
% figure
% subplot(3,1,1)
% plot(T,Om2(1,:),T,euler.Om2(1,:))
% title('\Omega_{2x}')
% subplot(3,1,2)
% plot(T,Om2(2,:),T,euler.Om2(2,:))
% title('\Omega_{2y}')
% subplot(3,1,3)
% plot(T,Om2(3,:),T,euler.Om2(3,:))
% title('\Omega_{2z}')
% legend('S2','Euler')
% 


% figure
% subplot(3,1,1)
% plot(T,Om1(1,:)-euler.Om1(1,:))
% title('\Omega_{1x} error')
% subplot(3,1,2)
% plot(T,Om1(2,:)-euler.Om1(2,:))
% title('\Omega_{1y} error')
% subplot(3,1,3)
% plot(T,Om1(3,:)-euler.Om1(3,:))
% title('\Omega_{1z} error')
% 
% figure
% subplot(3,1,1)
% plot(T,Om2(1,:)-euler.Om2(1,:))
% title('\Omega_{2x} error')
% subplot(3,1,2)
% plot(T,Om2(2,:)-euler.Om2(2,:))
% title('\Omega_{2y} error')
% subplot(3,1,3)
% plot(T,Om2(3,:)-euler.Om2(3,:))
% title('\Omega_{2z} error')
% 
% figure
% subplot(3,1,1)
% plot(T,KE - euler.KE)
% title('Kinetic Energy Error')
% subplot(3,1,2)
% plot(T,PE-euler.PE)
% title('Potential Energy Error')
% subplot(3,1,3)
% plot(T,TE-euler.TE)
% title('Total Energy Error')



function dx = acrobot_sim_s2(t,x,data)

m1 = data.m;
m2 = data.m;
l1 = data.l;
l2 = data.l;
g  = data.g;
e3 = data.e3;

q1 = x(1:3);
q2 = x(4:6);
Om1 = x(7:9);
Om2 = x(10:12);

dq1 = cross2(Om1,q1);
dq2 = cross2(Om2,q2);

J = [(m1+m2)*l1^2*eye(3) -(1/2)*m2*l1*l2*hat(q1)*hat(q2);
    -(1/2)*m2*l1*l2*hat(q2)*hat(q1) m2*l2^2*eye(3)];
C = [-(1/2)*m2*l1*l2*norm(Om2)^2*hat(q1)*q2;-(1/2)*m2*l1*l2*norm(Om1)^2*hat(q2)*q1];
G = [(m1+m2)*g*l1*hat(q1)*e3;m2*g*l2*hat(q2)*e3];
B =  [zeros(3);hat(q2)];

u = [0;0;0];

dOm = J\(B*u-G-C); 

dx =[dq1;dq2;dOm];

end


function dx = acrobot_sim_euler(t,x,data)

m1 = data.m;
m2 = data.m;
mt = data.mt;
l1 = data.l;
l2 = data.l;
g  = data.g;

ph1 = x(1);
th1 = x(2);
ph2 = x(3);
th2 = x(4);
dph1 = x(5);
dth1 = x(6);
dph2 = x(7);
dth2 = x(8);

D = acrobot_inertia_matrix(l1,l2,m1,m2,ph1,ph2,th1,th2);
H = acrobot_coriolis_and_gravity(dph1,dph2,dth1,dth2,g,l1,l2,m1,m2,ph1,ph2,th1,th2);
dx(1:4) = [dph1;dth1;dph2;dth2];

B = [-1 0;0 -1;1 0;0 1];
u = [0;0];

out = D\(H+B*u);
dx(5:8)=out;
dx=dx';
end 

function dx = acrobot_sim_so3(t,x,data)

m1 = data.m;
m2 = data.m;

l1 = data.l*data.e3;
lc1 = data.l*data.e3;
l2 = data.l*data.e3;
lc2 = data.l*data.e3;

J1 = data.J;
J2 = data.J;

g = data.g;
e3 = data.e3;

R1w = reshape(x(1:9),3,3);
omega1w = x(10:12);
R2w = reshape(x(13:21),3,3);
omega2w = x(22:24);

% converting R's and omegas from world-frame to base-frame
R1 = R1w';
R2 = R2w';
omega1 = R1*omega1w; % Here R1 and R2 are defined w.r.t world-frame
omega2 = R2*omega2w; 

D = [(J1-m1*hat(lc1)^2-m2*hat(l1)^2) -m2*hat(l1)*R1'*R2*hat(lc2);
     -m2*hat(lc2)*R2'*R1*hat(l1) (J2-m2*hat(lc2)^2)];

C = [hat(omega1)*(J1-(m1*hat(lc1)^2)-(m2*hat(l1)^2))*omega1 + m2*hat(l1)*R1'*R2*hat(omega2)^2*lc2;
     hat(omega2)*(J2-(m2*hat(lc2)^2))*omega2 + m2*hat(lc2)*R2'*R1*hat(omega1)^2*l1];

G = [(m1*g*hat(lc1)*R1'*e3 + m2*g*hat(l1)*R1'*e3) ;
     (m2*g*hat(lc2)*R2'*e3)];

B = [-(R1'*R2);eye(3)];  

M = [0;0;0];

out = D\(-C -G + B*M);

% R1_dot = R1*hat(omega1);
% R2_dot = R2*hat(omega2);
% 
% omega1_dot = out(1:3);
% omega2_dot = out(4:6);

% converting R's and omegas from base-frame to world-frame
R1_dot = -hat(omega1)*R1w;
R2_dot = -hat(omega2)*R2w;

omega1_dot = R1*out(1:3);
omega2_dot = R2*out(4:6);

dx =[reshape(R1_dot,9,1);omega1_dot;reshape(R2_dot,9,1);omega2_dot];

end