%function []= Crank_Slider_Analyzer()

%Created by Jack Manifold
%Created on March 7, 2021

%% Inputs

%Testing values
O_2 = [0,0];
a = 4;                      %Length of input
b = 4;                      %Length of coupler
c = 2.59;                   %Length of offset
omega_2 = 3;                %rads/sec
alpha_2 = 3;                %rad/sec^2
s1 = 15;                    %distance from first point of interest from O2
delta_2 = -19.3*pi/180;     %radians
s2 = 12;                    %distance from second point of interest from O2
s3 = 7;                     %distance of intermediate length from O2

%Setting up the range of Theta 2
theta2_min = 44.89*pi/180;       %Minimum angle of theta 2
theta2_max = 121.45*pi/180;      %Maximum angle of theta 2

theta2_increment = 0.05;      %The stepping factor of theta 2

%% Inital calculation and Preallocating
theta_2 = theta2_min:theta2_increment:theta2_max;


%Setting up matrixes
theta_3 = zeros(length(theta_2),2);
d = zeros(length(theta_2),2);
omega_3 = zeros(length(theta_2),2);

d_dot = zeros(length(theta_2),2);

Va_x = zeros(length(theta_2),2);
Va_y = zeros(length(theta_2),2);

Vba_x = zeros(length(theta_2),2);
Vba_y = zeros(length(theta_2),2);

Vb_x = zeros(length(theta_2),2);
Vb_y = zeros(length(theta_2),2);

Va = zeros(length(theta_2),2);
Vba = zeros(length(theta_2),2);
Vb = zeros(length(theta_2),2);

Vs1_x = zeros(length(theta_2),1);
Vs1_y = zeros(length(theta_2),1);
Vs1 = zeros(length(theta_2),1);

Vs2_x = zeros(length(theta_2),1);
Vs2_y = zeros(length(theta_2),1);
Vs2 = zeros(length(theta_2),1);

alpha_3 = zeros(length(theta_2),2);
d_2dot = zeros(length(theta_2),2);

Aa = zeros(length(theta_2),3);
Ab = zeros(length(theta_2),2);

As_1 = zeros(length(theta_2),3);
As_2 = zeros(length(theta_2),3);

Pt_A = zeros(length(theta_2),2);
Pt_B = zeros(length(theta_2),4);

Ps_1 = zeros(length(theta_2),2);
Ps_2 = zeros(length(theta_2),2);
Ps_3 = zeros(length(theta_2),2);

%% Calculating Values for Positions, Velocities, and Accelerations
for i = 1:length(theta_2)
    
    %Position Calculations
    
    theta_3(i,1) = asin((a*sin(theta_2(i))-c)/b);  %angle < pi/2
    theta_3(i,2) = asin(-(a*sin(theta_2(i))-c)/b)+pi; %angle > pi/2
    
    %Length of d or the distance of the crusher block from O2
    d(i,1) = a*cos(theta_2(i))-b*cos(theta_3(i,1)); %angle < pi/2
    d(i,2) = a*cos(theta_2(i))-b*cos(theta_3(i,2)); %angle > pi/2
    
    %Positions of Pt A and Pt B
    Pt_A(i,1) = a*cos(theta_2(i));
    Pt_A(i,2) = a*sin(theta_2(i));
    
    Pt_B(i,1) = a*cos(theta_2(i)) + b*cos(theta_3(i,1)); %X angle > pi/2
    Pt_B(i,2) = c;                                       %Y
    Pt_B(i,3) = a*cos(theta_2(i)) + b*cos(theta_3(i,2)); %X angle < pi/2
    Pt_B(i,4) = c;                                       %Y 
    
    %Position of Point S1
    Ps_1(i,1) = s1 * cos(theta_2(i)+delta_2);
    Ps_1(i,2) = s1 * sin(theta_2(i)+delta_2);
    
    %Position of Point S2
    Ps_2(i,1) = s2 * cos(theta_2(i));
    Ps_2(i,2) = s2 * sin(theta_2(i));
    
    %Intermediate position used in plotting
    Ps_3(i,1) = s3 * cos(theta_2(i));
    Ps_3(i,2) = s3 * sin(theta_2(i));
    
    
    %Velocity Calculations
    
    %Angular Velocity for Omega 3
    omega_3(i,1) = (a*omega_2*cos(theta_2(i)))/(b*cos(theta_3(i,1))); %angle < pi/2
    omega_3(i,2) = (a*omega_2*cos(theta_2(i)))/(b*cos(theta_3(i,2))); %angle > pi/2
    
    %Speed of point b in relation to O2
    d_dot(i,1) = -a*omega_2*sin(theta_2(i)) + b*omega_3(i,1)*sin(theta_3(i,1)); %angle < pi/2
    d_dot(i,2) = -a*omega_2*sin(theta_2(i)) + b*omega_3(i,2)*sin(theta_3(i,2)); %angle > pi/2
           
    %Velocity of Point A
    Va_x(i,1) = a*omega_2*(-sin(theta_2(i)));  %angle < pi/2
    Va_y(i,1) = a*omega_2*cos(theta_2(i));     
    
    Va_x(i,2) = a*omega_2*(-sin(theta_2(i)));  %angle > pi/2
    Va_y(i,2) = a*omega_2*cos(theta_2(i)); 
    
    Va(i,1) = ((Va_x(i,1))^2 + (Va_y(i,1))^2)^0.5;  %angle < pi/2
    Va(i,2) = ((Va_x(i,2))^2 + (Va_y(i,2))^2)^0.5;  %angle > pi/2
    
    %Velocity of Point B relative to A
    Vba_x(i,1) = -b*omega_3(i,1)*(-sin(theta_3(i,1))); %angle < pi/2
    Vba_y(i,1) = -b*omega_3(i,1)*(cos(theta_3(i,1)));
    
    Vba_x(i,2) = -b*omega_3(i,2)*(-sin(theta_3(i,2))); %angle > pi/2
    Vba_y(i,2) = -b*omega_3(i,2)*(cos(theta_3(i,2)));
    
    Vba(i,1) = ((Vba_x(i,1))^2 + (Vba_y(i,1))^2)^0.5; %angle < pi/2
    Vba(i,2) = ((Vba_x(i,2))^2 + (Vba_y(i,2))^2)^0.5; %angle > pi/2
    
    %Velocity of Point B
    Vb_x(i,1) = Va_x(i,1)+Vba_x(i,1);   %angle < pi/2
    Vb_y(i,1) = Va_y(i,1)+Vba_y(i,1);
    
    Vb_x(i,2) = Va_x(i,2)+Vba_x(i,2);   %angle > pi/2
    Vb_y(i,2) = Va_y(i,2)+Vba_y(i,2);
    
    Vb(i,1) = ((Vb_x(i,1))^2+(Vb_y(i,1))^2)^0.5;  %angle < pi/2
    Vb(i,2) = ((Vb_x(i,2))^2+(Vb_y(i,2))^2)^0.5;  %angle > pi/2
    
    %Velocity of Points Vs1 and Vs2
    
    Vs1_x(i) = -s1*omega_2*sin(theta_2(i)+delta_2);
    Vs1_y(i) =  s1*omega_2*cos(theta_2(i)+delta_2);
    
    Vs1(i) = ((Vs1_x(i))^2+(Vs1_y(i))^2)^0.5;
    
    Vs2_x(i) = -s2*omega_2*sin(theta_2(i));
    Vs2_y(i) =  s2*omega_2*cos(theta_2(i));
    
    Vs2(i) = ((Vs2_x(i))^2+(Vs2_y(i))^2)^0.5;
    
    
    %Acceleration Calculations
    
    
    %Angular Acceleration for Alpha 3
        %when theta 3 < pi/2
    alpha_3(i,1) = (a*alpha_2*cos(theta_2(i))...
        -a*(omega_2)^2*sin(theta_2(i))...
        +b*(omega_3(i,1))^2*sin(theta_3(i,1)))/(b*cos(theta_3(i,1)));
    
        %when theta 3 > pi/2
    alpha_3(i,2) = (a*alpha_2*cos(theta_2(i))...
        -a*(omega_2)^2*sin(theta_2(i))...
        +b*(omega_3(i,2))^2*sin(theta_3(i,2)))/(b*cos(theta_3(i,2)));
    
    %Acceleration of Length d
        %when theta 3 < pi/2
    d_2dot(i,1) = (b.*alpha_3(i,1)*sin(theta_3(i,1))...
        +b*(omega_3(i,1))^2*cos(theta_3(i,1)))...
        -(a*alpha_2*sin(theta_2(i))+a*(omega_2)^2*cos(theta_2(i)));
    
        %when theta 3 > pi/2
    d_2dot(i,2) = (b*alpha_3(i,2)*sin(theta_3(i,2))...
        +b*(omega_3(i,2))^2*cos(theta_3(i,2)))...
        -(a*alpha_2*sin(theta_2(i))+a*(omega_2)^2*cos(theta_2(i)));
    
    %Acceleration of Points A and B
    Aa(i,1) = -a*alpha_2*sin(theta_2(i))-a*(omega_2)^2*cos(theta_2(i)); %normal vector
    Aa(i,2) =  a*alpha_2*cos(theta_2(i))-a*(omega_2)^2*sin(theta_2(i)); %tangential vector
    Aa(i,3) = Aa(i,1)+Aa(i,2);
    
    Ab(i,1) = d_2dot(i,1); %tangential acceleration when theta 3 < pi/2
    Ab(i,2) = d_2dot(i,2); %tangential acceleration when theta 3 > pi/2
    
    %Ab sees no normal acceleration!
    
    %Acceleration of Points of interest
    As_1(i,1) = -s1 * alpha_2*sin(theta_2(i)+delta_2)...
        -s1*(omega_2)^2*cos(theta_2(i)+delta_2);
    
    As_1(i,2) =  s1 * alpha_2*cos(theta_2(i)+delta_2)...
        -s1*(omega_2)^2*sin(theta_2(i)+delta_2);
    
    As_1(i,3) = As_1(i,1)+As_1(i,2); %magnitude of acceleration for POI 1
    
    
    As_2(i,1) = -s2 * alpha_2*sin(theta_2(i))...
        -s2*(omega_2)^2*cos(theta_2(i)); %tangential acceleration
    
    As_2(i,2) =  s2 * alpha_2*cos(theta_2(i))...
        -s2*(omega_2)^2*sin(theta_2(i)); %normal acceleration
    
    As_2(i,3) = As_2(i,1)+As_2(i,2); %magnitude of acceleration for POI 2
    
end

%% Plotting

figure()
set(gcf,'position',[100 100 500 500])
scatter(O_2(1),O_2(2),140,'b')
hold on

%Plot Original Position
scatter(Pt_A(1,1),Pt_A(1,2),140,'b')
hold on
scatter(Pt_B(1,1),Pt_B(1,2),140,'b')
hold on
scatter(Ps_1(1,1),Ps_1(1,2),140,'b')
hold on
scatter(Ps_2(1,1),Ps_2(1,2),140,'b')
hold on
text([Pt_A(1,1)-2.3,Pt_B(1,1),Ps_1(1,1),Ps_2(1,1)-1],[Pt_A(1,2),Pt_B(2,2)-1,Ps_1(2,2)+1,Ps_2(2,2)+1],{'Pt A','Pt B','POI 1','POI 2'})

line([O_2(1),Pt_A(1,1)],[O_2(2),Pt_A(1,2),],'linewidth',5)
line([Pt_A(1,1),Pt_B(1,1)],[Pt_A(1,2),Pt_B(1,2)],'linewidth',5)
line([Pt_A(1,1),Ps_2(1,1)],[Pt_A(1,2),Ps_2(1,2)],'linewidth',5,'Linestyle','--')
line([Pt_A(1,1),Ps_3(1,1)],[Pt_A(1,2),Ps_3(1,2)],'linewidth',5,'Linestyle','-.')
line([Ps_3(1,1),Ps_1(1,1)],[Ps_3(1,2),Ps_1(1,2)],'linewidth',5,'Linestyle','-.')
%hold on 

plot(Pt_A(:,1),Pt_A(:,2),':')
plot(Pt_B(:,1),Pt_B(:,2),':')
plot(Ps_1(:,1),Ps_1(:,2),':')
plot(Ps_2(:,1),Ps_2(:,2),':')

xlim([-7,15])
ylim([-5,20])

hold on
title('Position Plot for Mechanism''s Range of Motion')
xlabel('X positions in inches')
ylabel('Y positions in inches')


figure();
set(gcf,'Position',[100 100 500 500])
for i = 1:length(theta_2)
    scatter(O_2(1),O_2(2),100,'b')
    hold on
    
    scatter(Pt_A(i,1),Pt_A(i,2),100,'b')
    hold on
    scatter(Pt_B(i,1),Pt_B(i,2),100,'b')
    hold on
    scatter(Ps_1(i,1),Ps_1(i,2),100,'b')
    hold on 
    scatter(Ps_2(i,1),Ps_2(i,2),100,'b')
    hold on
    
    
    line([O_2(1),Pt_A(i,1)],[O_2(2),Pt_A(i,2),],'linewidth',5)
    line([Pt_A(i,1),Pt_B(i,1)],[Pt_A(i,2),Pt_B(i,2)],'linewidth',5)
    line([Pt_A(i,1),Ps_3(i,1)],[Pt_A(i,2),Ps_3(i,2)],'linewidth',5,'Linestyle','-.')
    line([Ps_3(i,1),Ps_1(i,1)],[Ps_3(i,2),Ps_1(i,2)],'linewidth',5,'Linestyle','-.')
    line([Pt_A(i,1),Ps_2(i,1)],[Pt_A(i,2),Ps_2(i,2)],'linewidth',5,'Linestyle','-.')
    
    plot(Pt_A(:,1),Pt_A(:,2),':')
    plot(Pt_B(:,1),Pt_B(:,2),':')
    plot(Ps_1(:,1),Ps_1(:,2),':')
    plot(Ps_2(:,1),Ps_2(:,2),':')

    xlim([-7,15])
    ylim([-5,20])
   
    title('Position Plot for Mechanism''s Range of Motion')
    xlabel('X positions in inches')
    ylabel('Y positions in inches')

    hold off
    M(i) = getframe;
    
    if i == length(theta_2)
        t = i;
        count = i;
        while t ~= 0
            scatter(O_2(1),O_2(2),100,'b')
            hold on

            scatter(Pt_A(t,1),Pt_A(t,2),100,'b')
            hold on
            scatter(Pt_B(t,1),Pt_B(t,2),100,'b')
            hold on
            scatter(Ps_1(t,1),Ps_1(t,2),100,'b')
            hold on
            scatter(Ps_2(t,1),Ps_2(t,2),100,'b')
            hold on
            
            line([O_2(1),Pt_A(t,1)],[O_2(2),Pt_A(t,2),],'linewidth',5)
            line([Pt_A(t,1),Pt_B(t,1)],[Pt_A(t,2),Pt_B(t,2)],'linewidth',5)
            line([Pt_A(t,1),Ps_2(t,1)],[Pt_A(t,2),Ps_2(t,2)],'linewidth',5,'Linestyle','-.')
            line([Pt_A(t,1),Ps_3(t,1)],[Pt_A(t,2),Ps_3(t,2)],'linewidth',5,'Linestyle','-.')
            line([Ps_3(t,1),Ps_1(t,1)],[Ps_3(t,2),Ps_1(t,2)],'linewidth',5,'Linestyle','-.')

            plot(Pt_A(:,1),Pt_A(:,2),':')
            plot(Pt_B(:,1),Pt_B(:,2),':')
            plot(Ps_1(:,1),Ps_1(:,2),':')
            plot(Ps_2(:,1),Ps_2(:,2),':')

            xlim([-7,15])
            ylim([-5,20])

            title('Position Plot for Mechanism''s Range of Motion')
            xlabel('X positions in inches')
            ylabel('Y positions in inches')
            
            hold off
            M(count) = getframe;

            t = t-1;
            
            count = count +1;
        end
    end
    
end

movie(M,10,60)







