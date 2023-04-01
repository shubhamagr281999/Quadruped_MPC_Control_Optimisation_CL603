clc
clear
%% defining trajectory
S=[0,3/4,2/4,1/4];  % defines the gait pattern
dutyFactor=1/4;
Cycle_time=2;
stance_time=Cycle_time*dutyFactor;
swing_time=Cycle_time*(1-dutyFactor);
Lstance=0.1;
Lswing=Lstance;
nominal_height=0.085;
%% swing trjacetory
Xs=-Lswing/2:Lswing/((1-dutyFactor)*points_per_Cycle-1):Lswing/2+0.004;
Zs=sin((Xs+Lswing/2)*pi/Lswing)-nominal_height;
Ys=0.05;
%% stance trjectory
Xst=-Lstance/2:Lstance/(dutyFactor*points_per_Cycle-1):Lstance/2+0.004;
Yst=0.05;
Zst=nominal_height;
%% robot_setup for inverse kinematics



