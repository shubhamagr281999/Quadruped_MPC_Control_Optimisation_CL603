%% assignment of DH paramaters to be used by the function built and in link formation
l1=0.1;
l2=0.2092;
l3=0.3294;
l4=0.2;
alpha=[pi/2,0,0,0];
d=[0,0,0,0];
a=[l1,l2,l3,l4];
tor=zeros(21,21,4);

%% defining link and initializing bot 
L1=Link('revolute','d',d(1,1),'a',a(1,1),'alpha',alpha(1,1));
L2=Link('revolute','d',d(1,2),'a',a(1,2),'alpha',alpha(1,2));
L3=Link('revolute','d',d(1,3),'a',a(1,3),'alpha',alpha(1,3));
L4=Link('revolute','d',d(1,4),'a',a(1,4),'alpha',alpha(1,4));
bot=SerialLink([L1,L2,L3,L4],'name','hexapod');
%% comparing the forward dynamincs with our code and fkine
% T=ForwardKinematics([pi,pi/2,pi/4,pi/3,pi/2,pi/6])
theta=[0,1.1,-2.1,0];
theta(4)=-pi/2-theta(2)-theta(3);
T=bot.fkine(theta)
%% code for torque
for i=1:1:31
    theta=[i*0.1,pi/3,-pi*0.8,-pi/2-pi/3+pi*0.8];
    Jacob=bot.jacob0(theta);
    tor(i,:)=(Jacob'*[0;78.4;196;0;0;0])'
end
%% plots
figure
plot(0.1:0.1:3.1,tor(:,2))

%% theta2,3 vary
for i=1:1:14
    for j=1:1:21
        theta=[0.5,i*0.1,j*-0.1,-pi/2-0.1*i+0.1*j];
        Jacob=bot.jacob0(theta);
        tor(i,j,:)=(Jacob'*[0;78.4;196;0;0;0])';
    end
end
%% plot2
figure
surf(0.1:0.1:2.1,0.1:0.1:2.1,abs(tor(:,:,2)))
xlabel('A2')
ylabel('A3')
%% inverse kinematics
eff_pose=[0.4242,0.2449,0.2375];
t0_3=[eff_pose(1),eff_pose(2),eff_pose(3)+l4];
theta=zeros(4,1);
theta(1)=atan2(t0_3(2),t0_3(1));

t0_3(1)=t0_3(1)-l1*cos(theta(1));
t0_3(2)=t0_3(2)-l1*sin(theta(1));

x = optimvar('x',2);
eq1=sqrt(t0_3(1).^2+t0_3(2).^2)-l2*cos(x(1))-l3*cos(x(1)+x(2))==0;
eq2=t0_3(3)-l2*sin(x(1))-l3*sin(x(1)+x(2))==0
prob=eqnproblem;
prob.Equations.eq1 = eq1;
prob.Equations.eq2 = eq2;
x0.x = [0.1 0.05];
[sol,fval,exitflag] = solve(prob,x0);
theta(2)=sol.x(1);
theta(3)=sol.x(2);
theta(4)=-pi/2-theta(2)-theta(3);