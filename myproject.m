tic;
%%
% defining constants and variable matrices:
%
thetta_3 = [];
thetta_4 = [];
r_5 = [];
r_6 = [];
y = [];
% consants :
r_1 = 720.00/1000.0;
r_2 = 130.00/1000.0;  
r_3 = 640.00/1000.0;
r_4 = 260.00/1000.0;
r_7 = 550.00/1000.0;
r_8 = 400.00/1000.0;

thetta_1 = 0.0;%defined in degrees
thetta_8 = 0.0;
thetta_53 = 145.0;
thetta_65 = 330.0;
thetta_73 = 120.0;

% inputs:
omega_2 = 2.0*pi;%rad per sec
alpha_2 = 0;%rad per square sec

% Set option to display information after each iterationc
options=optimset('Display','iter');
% answers :
firstLoopPositionAnswer = [];
firstLoopVelocityAnswer = [];
firstLoopAccelerationAnswer = [];
secondLoopPositionAnswer = [];
secondLoopVelocityAnswer = [];
secondLoopAccelerationAnswer = [];
%%
% position solution
%

%
% First loop first
%
% thetta_3 = x(1), thetta_4 = x(2), thetta_2=k 

for k = 0:360
firstLoopPositionSolution = @(x) [ r_2*cosd(k) + r_3*cosd(x(1)) - r_4*cosd(x(2)) - r_1*cosd(thetta_1);...
                                   r_2*sind(k) + r_3*sind(x(1)) - r_4*sind(x(2)) - r_1*sind(thetta_1)];
% Make a starting guess at the solution
x0=[30;30];

% Solve the system
[x,fval,exitflag] = fsolve(firstLoopPositionSolution, x0, options);

firstLoopPositionAnswer = cat(2, firstLoopPositionAnswer, x);
end
thetta_3 = firstLoopPositionAnswer(1, :);
thetta_4 = firstLoopPositionAnswer(2, :);
thetta_5 = thetta_3 + thetta_53;
thetta_6 = thetta_5 + thetta_65;
thetta_7 = thetta_3 + thetta_73;
clc
%
% Second loop:
% y(1) = r_5, y(2) = r_6
for k = 1:361
y = ([cosd(thetta_5(k)) cosd(thetta_6(k));sind(thetta_5(k)) sind(thetta_6(k))])\...
    [+r_7*cosd(thetta_7(k)) - r_8*cosd(thetta_8) + r_4*cosd(thetta_4(k));...
     +r_7*sind(thetta_7(k)) - r_8*sind(thetta_8) + r_4*sind(thetta_4(k))];
secondLoopPositionAnswer = cat(2, secondLoopPositionAnswer, y);
end
r_6 = secondLoopPositionAnswer(2, :);
r_5 = secondLoopPositionAnswer(1, :);

%%
% Velocity Solution:
%

%
% First loop first
%
% omega_3 = x(1), omega_4 = x(2), thetta_2=k 
x=[];
for k = 1:361
firstLoopVelocitySolution = @(x) [ -r_2*omega_2*sind(k-1) - r_3*x(1)*sind(thetta_3(k)) + r_4*x(2)*sind(thetta_4(k));...
                                   +r_2*omega_2*cosd(k-1) + r_3*x(1)*cosd(thetta_3(k)) - r_4*x(2)*cosd(thetta_4(k))];
% Make a starting guess at the solution
x0=[20;20];

% Solve the system
[x,fval,exitflag] = fsolve(firstLoopVelocitySolution, x0, options);

firstLoopVelocityAnswer = cat(2, firstLoopVelocityAnswer, x);
end
omega_3 = firstLoopVelocityAnswer(1, :);
omega_4 = firstLoopVelocityAnswer(2, :);
omega_5 = omega_3;
omega_6 = omega_3;
omega_7 = omega_3;
y= [];
%
% Second loop:
% y(1) = v_5, y(2) = v_6
for k = 1:361
y = ([-cosd(thetta_5(k)) -cosd(thetta_6(k));-sind(thetta_5(k)) -sind(thetta_6(k))])\...
    [r_4*omega_4(k)*sind(thetta_4(k)) + r_7*omega_7(k)*sind(thetta_7(k))...
     - r_5(k)*omega_5(k)*sind(thetta_5(k))- r_5(k)*omega_5(k)*sind(thetta_5(k))...
     - r_6(k)*omega_6(k)*sind(thetta_6(k))- r_6(k)*omega_6(k)*sind(thetta_6(k));...
     -r_4*omega_4(k)*cosd(thetta_4(k)) - r_7*omega_7(k)*cosd(thetta_7(k))...
    + r_5(k)*omega_5(k)*cosd(thetta_5(k)) + r_5(k)*omega_5(k)*cosd(thetta_5(k))...
    + r_6(k)*omega_6(k)*cosd(thetta_6(k)) + r_6(k)*omega_6(k)*cosd(thetta_6(k))];
secondLoopVelocityAnswer = cat(2, secondLoopVelocityAnswer, y);
end
v_6 = secondLoopVelocityAnswer(2, :);
v_5 = secondLoopVelocityAnswer(1, :);

%%
% Acceleration Solution:
%

%
% First loop first
%
% alpha_3 = x(1), alpha_4 = x(2), thetta_2=k-1 
x=[];
for k = 1:361
firstLoopAccelerationSolution = @(x) [ -r_2*alpha_2*sind(k-1) - r_2*(omega_2^2)*cosd(k-1) ...
                                       - r_3*x(1)*sind(thetta_3(k)) - r_3*(omega_3(k)^2)*cosd(thetta_3(k))...
                                       + r_4*x(2)*sind(thetta_4(k)) + r_4*(omega_4(k)^2)*cosd(thetta_4(k));...
                                       % ? first eq. 
                                       +r_2*alpha_2*cosd(k-1) - r_2*(omega_2^2)*sind(k-1) ...
                                       + r_3*x(1)*cosd(thetta_3(k)) - r_3*(omega_3(k)^2)*sind(thetta_3(k))...
                                       - r_4*x(2)*cosd(thetta_4(k)) + r_4*(omega_4(k)^2)*sind(thetta_4(k));...
                                       % ? second eq.
                                     ];
% Make a starting guess at the solution
x0=[1;1];

% Solve the system
[x,fval,exitflag] = fsolve(firstLoopAccelerationSolution, x0, options);

firstLoopAccelerationAnswer = cat(2, firstLoopAccelerationAnswer, x);
end
alpha_3 = firstLoopAccelerationAnswer(1,:);
alpha_4 = firstLoopAccelerationAnswer(2,:);
alpha_5 = alpha_3;
alpha_6 = alpha_3;
alpha_7 = alpha_3;
y=[];
%
% Second loop:
% y(1) = a_5, y(2) = a_6
for k = 1:361
y = ([-cosd(thetta_5(k)) -cosd(thetta_6(k));-sind(thetta_5(k)) -sind(thetta_6(k))])\...
    [r_4*alpha_4(k)*sind(thetta_4(k)) + r_4*(omega_4(k)^2)*cosd(thetta_4(k)) + r_7*alpha_7(k)*sind(thetta_7(k)) + r_7*(omega_7(k)^2)*cosd(thetta_7(k)) - 2.0*v_5(k)*omega_5(k)*sind(thetta_5(k)) - r_5(k)*alpha_5(k)*sind(thetta_5(k)) + r_5(k)*(omega_5(k)^2)*cosd(thetta_5(k)) - 2.0*v_6(k)*omega_6(k)*sind(thetta_6(k)) - r_6(k)*alpha_6(k)*sind(thetta_6(k)) + r_6(k)*(omega_6(k)^2)*cosd(thetta_6(k));...
     -r_4*alpha_4(k)*cosd(thetta_4(k)) + r_4*(omega_4(k)^2)*sind(thetta_4(k)) - r_7*alpha_7(k)*cosd(thetta_7(k)) + r_7*(omega_7(k)^2)*sind(thetta_7(k)) + 2.0*v_5(k)*omega_5(k)*cosd(thetta_5(k)) + r_5(k)*alpha_5(k)*cosd(thetta_5(k)) - r_5(k)*(omega_5(k)^2)*cosd(thetta_5(k)) + 2.0*v_6(k)*omega_6(k)*cosd(thetta_6(k)) + r_6(k)*alpha_6(k)*cosd(thetta_6(k)) - r_6(k)*(omega_6(k)^2)*cosd(thetta_6(k))];
secondLoopAccelerationAnswer = cat(2, secondLoopAccelerationAnswer, y);
end
a_6 = secondLoopAccelerationAnswer(2, :);
a_5 = secondLoopAccelerationAnswer(1, :);

%%
% plotting results:
%
plot(0:360, thetta_6-360);
hold on;
xlabel('\theta_{2}');
ylabel('\theta_{6}');
title('Position');
grid on;
saveas(gcf, 'position', 'bmp');
saveas(gcf, 'position', 'fig');
close(gcf);

%velocity
plot(0:360, omega_6);
hold on;
xlabel('\theta_{2}');
ylabel('\omega_{6}');
title('Velocity');
grid on;
saveas(gcf, 'velocity', 'bmp');
saveas(gcf, 'velocity', 'fig');
close(gcf);

%accerelation
plot(0:360, alpha_6);
hold on;
xlabel('\theta_{2}');
ylabel('\alpha_{6}');
title('Accerelation');
grid on;
saveas(gcf, 'accerelation', 'bmp');
saveas(gcf, 'accerelation', 'fig');
close(gcf);
clc;
fprintf("All Done!\nFinished in "+toc+" seconds\n");
tic;
%% ************************************************
% Calculating Kintics and Inertial forces:
%
% AX=B

a_mat = zeros(15);
b_mat = zeros(15, 1);
f_12_x = [];
f_12_y = [];
f_23_x = [];
f_23_y = [];
f_34_x = [];
f_34_y = [];
f_35_x = [];
f_35_y = [];
f_14_x = [];
f_14_y = [];
f_16_x = [];
f_16_y = [];
f_56_x = [];
f_56_y = [];
t_12 = [];

m_2 = 1.50;%kg
m_3 = 1.00;
m_4 = 0.7;
m_5 = 7.600;
m_6 = 2.6;

i_2 = 8.00;
i_3 = 10.00;
i_4 = 5.00;
i_5 = 6.600;
i_6 = 3.3;
% *************** 2
a_mat(1, 1) = 1.00;
a_mat(1, 3) = 1.00;
a_mat(2, 2) = 1.00;
a_mat(2, 4) = 1.00;
% 
% 3rd row of a_mat
% must be calculated in loop
%
a_mat(3, 15) = 1.00;
% *************** /2
% *************** 3 
a_mat(4, 3) = 1.00;
a_mat(4, 5) = -1.00;
a_mat(4, 7) = -1.00;

a_mat(5, 4) = 1.00;
a_mat(5, 6) = -1.00;
a_mat(5, 8) = 1.00;
% 
% 6th row of a_mat
% must be calculated in loop
%
% *************** /3
% *************** 4 
a_mat(7, 5) = -1.00;
a_mat(7, 9) = 1.00;

a_mat(8, 6) = -1.00;
a_mat(8, 10) = 1.00;
% 
% 9th row of a_mat
% must be calculated in loop
%
% *************** /4
% *************** 5 
a_mat(10, 7) = -1.00;
a_mat(10, 13) = -1.00;

a_mat(11, 8) = 1.00;
a_mat(11, 14) = 1.00;
% 
% 12th row of a_mat
% must be calculated in loop
%
% *************** /5
% *************** 6 
a_mat(13, 11) = 1.00;
a_mat(13, 13) = -1.00;

a_mat(14, 12) = 1.00;
a_mat(14, 14) = -1.00;
% 
% 15th row of a_mat
% must be calculated in loop
%
% *************** /6
r_12_x_mat = [];
r_12_y_mat = [];

r_32_x_mat = [];
r_32_y_mat = [];
co_r_3 = sqrt(3.00)/3.00;
R_5 = max (r_5);
R_6 = max (r_6);

for i = 1:361
    r_12_y = r_2/2.0*sind(i-1);
    r_12_x = r_2/2.0*cosd(i-1);
    r_32_y = r_2*sind(i-1)-r_12_y;
    r_32_x = r_2*cosd(i-1)-r_12_x;
    
    r_23_y = co_r_3*r_3*sind(thetta_3(i)+30.00);
    r_23_x = co_r_3*r_3*cosd(thetta_3(i)+30.00);
    r_43_y = co_r_3*r_3*sind(thetta_3(i)+150.00);
    r_43_x = co_r_3*r_3*cosd(thetta_3(i)+150.00);
    r_53_y = co_r_3*r_3*sind(thetta_3(i)+270.00);
    r_53_x = co_r_3*r_3*cosd(thetta_3(i)+270.00);
    
    r_14_y = r_4/2.0*sind(thetta_4(i));
    r_14_x = r_4/2.0*cosd(thetta_4(i));
    r_34_y = r_4*sind(thetta_4(i))-r_14_y;
    r_34_x = r_4*cosd(thetta_4(i))-r_14_x;
    
    r_65_y = R_5/2.0*sind(thetta_5(i));
    r_65_x = R_5/2.0*cosd(thetta_5(i));
    r_35_y = r_5(i)*sind(thetta_5(i))-r_65_y;
    r_35_x = r_5(i)*cosd(thetta_5(i))-r_65_x;
    
    r_16_y = R_6/2.0*sind(thetta_6(i));
    r_16_x = R_6/2.0*cosd(thetta_6(i));
    r_56_y = r_6(i)*sind(thetta_6(i))-r_16_y;
    r_56_x = r_6(i)*cosd(thetta_6(i))-r_16_x;
    
    a_mat(3, 1) = r_12_y;
    a_mat(3, 2) = -r_12_x;
    a_mat(3, 3) = -r_32_y;
    a_mat(3, 4) = -r_32_x;
    
    a_mat(6, 3) = -r_23_y;
    a_mat(6, 4) = r_23_x;
    a_mat(6, 5) = -r_43_y;
    a_mat(6, 6) = -r_43_x;
    a_mat(6, 7) = r_53_y;
    a_mat(6, 8) = r_53_x;
    
    a_mat(9, 5) = r_34_y;
    a_mat(9, 6) = r_34_x;
    a_mat(9, 9) = r_14_y;
    a_mat(9, 10) = r_14_x;
    
    a_mat(12, 7) = -r_35_y;
    a_mat(12, 8) = r_35_x;
    a_mat(12, 13) = r_65_y;
    a_mat(12, 14) = -r_65_x;

    a_mat(15, 11) = r_16_y;
    a_mat(15, 12) = r_16_x;
    a_mat(15, 13) = r_56_y;
    a_mat(15, 14) = r_56_x;
    
    
    a_g2_x = sqrt(r_12_x^2+r_12_y^2) * (alpha_2 * cosd(i-1+90) + omega_2^2 * cosd(i-1));
    a_g2_y = sqrt(r_12_x^2+r_12_y^2) * (alpha_2 * sind(i-1+90) + omega_2^2 * sind(i-1));
    
    a_o3_x = r_2 * (alpha_2 * cosd(i-1+90) + omega_2^2 * cosd(i-1));
    a_o3_y = r_2 * (alpha_2 * sind(i-1+90) + omega_2^2 * sind(i-1));
    
    a_g3_x = a_o3_x + sqrt(r_23_x^2+r_23_y^2) * (alpha_3(i) * cosd(thetta_3(i)+90) + omega_3(i)^2 * cosd(thetta_3(i)));
    a_g3_y = a_o3_y + sqrt(r_23_x^2+r_23_y^2) * (alpha_3(i) * sind(thetta_3(i)+90) + omega_3(i)^2 * sind(thetta_3(i)));
    
    
    a_g4_x = sqrt(r_14_x^2+r_14_y^2) * (alpha_4(i) * cosd(thetta_4(i)+90) + omega_4(i)^2 * cosd(thetta_4(i)));
    a_g4_y = sqrt(r_14_x^2+r_14_y^2) * (alpha_4(i) * sind(thetta_4(i)+90) + omega_4(i)^2 * sind(thetta_4(i)));
    
    
    a_g6_x = sqrt(r_16_x^2+r_16_y^2) * (alpha_6(i) * cosd(thetta_6(i)+90) + omega_6(i)^2 * cosd(thetta_6(i)));
    a_g6_y = sqrt(r_16_x^2+r_16_y^2) * (alpha_6(i) * sind(thetta_6(i)+90) + omega_6(i)^2 * sind(thetta_6(i)));
    a_o5_x = r_6(i) * (alpha_6(i) * cosd(thetta_6(i)+90) + omega_6(i)^2 * cosd(thetta_6(i)));
    a_o5_y = r_6(i) * (alpha_6(i) * sind(thetta_6(i)+90) + omega_6(i)^2 * sind(thetta_6(i)));
    a_g5_x = a_o5_x + sqrt(r_65_x^2+r_65_y^2) * (alpha_5(i) * cosd(thetta_5(i)+90) + omega_5(i)^2 * cosd(thetta_5(i)));
    a_g5_y = a_o5_y + sqrt(r_65_x^2+r_65_y^2) * (alpha_5(i) * sind(thetta_5(i)+90) + omega_5(i)^2 * sind(thetta_5(i)));
    
    b_mat(1, 1) = m_2 * a_g2_x;
    b_mat(2, 1) = m_2 * a_g2_y;
    b_mat(3, 1) = i_2 * alpha_2;
    
    b_mat(4, 1) = m_3 * a_g3_x;
    b_mat(5, 1) = m_3 * a_g3_y;
    b_mat(6, 1) = i_3 * alpha_3(i);
    
    b_mat(7, 1) = m_4 * a_g4_x;
    b_mat(8, 1) = m_4 * a_g4_y;
    b_mat(9, 1) = i_4 * alpha_4(i);
    
    b_mat(10, 1) = m_5 * a_g5_x;
    b_mat(11, 1) = m_5 * a_g5_y;
    b_mat(12, 1) = i_5 * alpha_5(i);
    
    b_mat(13, 1) = m_6 * a_g6_x;
    b_mat(14, 1) = m_6 * a_g6_y;
    b_mat(15, 1) = i_6 * alpha_6(i);
    
    inertial_forces = a_mat \ b_mat;
    
    f_12_x = cat(2, f_12_x, inertial_forces(1, 1));
    f_12_y = cat(2, f_12_y, inertial_forces(2, 1));
    
    f_23_x = cat(2, f_23_x, inertial_forces(3, 1));
    f_23_y = cat(2, f_23_y, inertial_forces(4, 1));
    
    f_34_x = cat(2, f_34_x, inertial_forces(5, 1));
    f_34_y = cat(2, f_34_y, inertial_forces(6, 1));
    
    f_35_x = cat(2, f_35_x, inertial_forces(7, 1));
    f_35_y = cat(2, f_35_y, inertial_forces(8, 1));
    
    f_14_x = cat(2, f_14_x, inertial_forces(9, 1));
    f_14_y = cat(2, f_14_y, inertial_forces(10, 1));
    
    f_16_x = cat(2, f_16_x, inertial_forces(11, 1));
    f_16_y = cat(2, f_16_y, inertial_forces(12, 1));
    
    f_56_x = cat(2, f_56_x, inertial_forces(13, 1));
    f_56_y = cat(2, f_56_y, inertial_forces(14, 1));
    
    t_12 = cat(2, t_12, inertial_forces(15, 1));
    
    
    
end
fprintf("All Done!\nFinished in "+toc+" seconds\n");
