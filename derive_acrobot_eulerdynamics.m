close all
addpath('/home/exx/Avinash/My Toolbox');

syms m1 m2 l1 l2 real
syms th1 th2 dth1 dth2 ddth1 ddth2 ph1 ph2 dph1 dph2 ddph1 ddph2 real
syms g real

%%% basis vectors in the world frame
e1 = [1;0;0];
e2 = [0;1;0];
e3 = [0;0;1];

%%% generalized coordinates defined as going from world-frame to base-frame
eu1 = [ph1;th1];
deu1 = [dph1;dth1];
ddeu1 = [ddph1;ddth1];

eu2 = [ph2;th2];
deu2 = [dph2;dth2];
ddeu2 = [ddph2;ddth2];

eu = [eu1;eu2];
deu = [deu1;deu2];
ddeu = [ddeu1;ddeu2];

%%% position vectors
q1 = Rx(ph1)*Ry(th1)*e3;
q2 = Rx(ph2)*Ry(th2)*e3;


%%% CoM positions
p1 = (l1)*q1;
p2 = (l1)*q1 + (l2)*q2;

%%% CoM velocities
dq1 = jacobian(q1,eu1)*deu1;
dq2 = jacobian(q2,eu2)*deu2;

dp1 = jacobian(p1,eu1)*deu1;
dp2 = jacobian(p2,eu1)*deu1 + jacobian(p2,eu2)*deu2;

%%% CoM accelarations
ddq1 = jacobian(dq1,eu1)*deu1 + jacobian(dq1,deu1)*ddeu1;  
ddq2 = jacobian(dq2,eu2)*deu2 + jacobian(dq2,deu2)*ddeu2;

%%% Energies 

KE = simplify(0.5*m1*dot(dp1,dp1) + 0.5*m2*dot(dp2,dp2));
PE = simplify(m1*g*p1(3) + m2*g*p2(3));


[D,C,G] = EoM(KE,PE,eu,deu);

D = simplify(D);
H = simplify(-C*deu-G) ; %% Coriolis and Gravity terms lumped

matlabFunction(KE,'file','acrobot_kinetic_energy');
matlabFunction(PE,'file','acrobot_potential_energy');
matlabFunction(D,'file','acrobot_inertia_matrix');
matlabFunction(H,'file','acrobot_coriolis_and_gravity');

%%%%%% Findings a map from  u_e to u_q %%%%

syms u_th12 u_ph12 %% euler actuations at the joint between q1 and q2

B = [-1 0;0 -1;1 0;0 1];
ddeu = simplify(D\(H + B*[u_ph12;u_th12])); 
ddeu1 = ddeu(1:2);
ddeu2 = ddeu(3:4);

dOm1 = simplify(hat(q1)*(jacobian(dq1,eu1)*deu1 + jacobian(dq1,deu1)*ddeu1));
dOm2 = simplify(hat(q2)*(jacobian(dq2,eu2)*deu2 + jacobian(dq2,deu2)*ddeu2));
Om1 = simplify(hat(q1)*dq1);
Om2 = simplify(hat(q2)*dq2);

J = [((m1/4)+m2+mt)*l1^2*eye(3) -(1/2)*m2*l1*l2*hat(q1)*hat(q2);-(1/2)*m2*l1*l2*hat(q2)*hat(q1) (1/4)*m2*l2^2*eye(3)];
C = [-(1/2)*m2*l1*l2*norm(Om2)^2*hat(q1)*q2;-(1/2)*m2*l1*l2*norm(Om1)^2*hat(q2)*q1];
G = [((1/2)*m1+m2+mt)*g*l1*hat(q1)*e3;(1/2)*m2*g*l2*hat(q2)*e3];

U_big = J*[dOm1;dOm2] + C + G; %simplify();
matlabFunction(U_big,'file','acrobot_input')