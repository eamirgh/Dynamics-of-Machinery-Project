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
[x,fval,exitflag] = fsolve(firstLoopPositionSolution, x0, options)

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
y = inv([cosd(thetta_5(k)) cosd(thetta_6(k));sind(thetta_5(k)) sind(thetta_6(k))])*...
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
[x,fval,exitflag] = fsolve(firstLoopVelocitySolution, x0, options)

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
y = inv([-cosd(thetta_5(k)) -cosd(thetta_6(k));-sind(thetta_5(k)) -sind(thetta_6(k))])*...
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
[x,fval,exitflag] = fsolve(firstLoopAccelerationSolution, x0, options)

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
y = inv([-cosd(thetta_5(k)) -cosd(thetta_6(k));-sind(thetta_5(k)) -sind(thetta_6(k))])*...
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
fprintf("All Done!\nFinished in "+toc+" seconds\n")