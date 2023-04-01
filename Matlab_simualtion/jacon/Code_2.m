clear
clc
close all

%% Step 1: Generate Data
%Animate a point moving along a 3D parametric curve
t = linspace(0,10,300);
% x = 5*cos(t);
% y = 2*sin(t);
% z = t;
traj = zeros(3,length(t));

%% Step 2: Draw/Render Scenario
figh = figure;
for k=1:length(t)
    %Clear the figure to start with a blank slate
    clf
    
%     %Extract data at the current time step
%     t_k = t(k);
%     x_k = x(k);
%     y_k = y(k);
%     z_k = z(k);
%     
%     %Where is the current point?
%     plot3(x_k, y_k, z_k, 'go', 'LineWidth', 3, 'MarkerSize', 15)
%     
%     %Plot the entire curve
%     hold on
%     plot3(x, y, z, 'b-', 'LineWidth', 2);
%     
%     %Add plotting options
%     grid on
%     xlabel('x')
%     ylabel('y')
%     zlabel('z')
%     title(['t = ',num2str(t_k)])
%     view([30 35])
% %     view([30+20*t_k 35])      %show how viewpoint can be manipulated

x = [0 1 1 0 0]+ k/100; ;
y = [0 0 1 1 0]+ k/100;
z = [0 0 0 0 0];


% hold on
plot3(x,y,z)
patch(x,y,z,'m')

x = [0 1 1 0 0]+ k/100; 
y = [0 0 0 0 0]+ k/100;
z = [0 0 1 1 0];

patch(x,y,z,'m')

x = [0 0 0 0 0]+ k/100;
y = [0 1 1 0 0]+ k/100;
z = [0 0 1 1 0];

patch(x,y,z,'m')

x = [0 1 1 0 0]+ k/100;
y = [0 0 1 1 0]+ k/100;
z = [1 1 1 1 1];  %%top face

x_tp = x;
y_tp = y;
z_tp = z;

patch(x,y,z,'m')

x = [0 1 1 0 0]+ k/100;
y = [1 1 1 1 1]+ k/100;
z = [0 0 1 1 0];

patch(x,y,z,'m')

x = [1 1 1 1 1]+ k/100;
y = [0 1 1 0 0]+ k/100;
z = [0 0 1 1 0];

patch(x,y,z,'m')

xlabel('x')
ylabel('y')
zlabel('z')

x = [0.5 0.5+cos(k/10)]+ k/100;
y = [0.5 0.5+sin(k/10)]+ k/100;
z = [1 1];
hold on

traj(1,k) = x(2);
traj(2,k) = y(2);
traj(3,k) = z(2);

plot3(traj(1,1:k),traj(2,1:k),traj(3,1:k),'k')
line(x,y,z)

xlim([-1 6])
ylim([-1 6])
zlim([0 6])

view([30 30])
% view([30+30*k/100 35])

grid on
    
    %% Step 3: Take a Snapshot
    %Save the frame
    movieVector(k) = getframe;
%     movieVector(k) = getframe(figh, [10 10 520 400]);   %manually specify getframe region
    
    %% Step 4: Advance Time
    %Happens automatically if using a for loop
end

%% Step 5: Save Movie
%Create a VideoWriter object and set properties
myWriter = VideoWriter('curve');            %create an .avi file
% myWriter = VideoWriter('curve','MPEG-4');   %create an .mp4 file
myWriter.FrameRate = 20;

%Open the VideoWriter object, write the movie, and close the file
open(myWriter);
writeVideo(myWriter, movieVector);
close(myWriter);

disp('DONE!')