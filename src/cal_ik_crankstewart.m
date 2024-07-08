% code based on the paper titled "stewart platform with fixed rotary actuators: a low cost design study"
% Input parameters: 
%       
%       r_base             [m] Radius of the circumcircle of the base hexagon. 
%                           The points where the servo arm is connected to
%                           the servo axis lie on this circle.
%
%       r_platform             [m] Radius of the circumcircle of the platform
%                           hexagon. The points where the rod is attached
%                           to the platform lie on this circle.
%       
%       servo_arm_l [m]    Length of the servo arm. All servo arms  
%                           should have the same length.
%
%       rod_l      [m] Length of the rod connecting the servo arm to
%                           the platform.
%       
%       alpha_B         [rad]   Half the angle between two servos arm    
%                           attachment points. The two servo's arms must  
%                           lie in the same plane.
%
%       alpha_P         [rad]   Half the angle between two rod attachment
%                           points. They also have to be on the same edge
%                           of the platform hexagon.
%
%       trans           [m] The vector which describes the wanted 
%                           translational movement of the platform. 
%                           [x, y, z]
%
%       orient          [rad]   The vector which describes the wanted
%                           rotational movement of the platform according
%                           to the three euler angles [ phi, theta, psi].
%                           
function [servo_angle_del]=cal_ik_crankstewart(r_base,r_platform, servo_arm_l, rod_l, alpha_B, alpha_P, trans, orient)
trans= trans(:);
orient= orient(:);
%% Define the Geometry of the platform
% Beta is the angle between the plane in which the servo arm moves and the
% xz-plane of the base CS.
beta= [ pi+pi/2, pi/2,...
        2*pi/3+pi+pi/2, 2*pi/3+pi/2,...
        4*pi/3+pi+pi/2, 4*pi/3+pi/2];
    
    %gamma_B represents the direction where the servo arm is attached to
    %the servo axis, we use polar coord system
    gamma_B= [alpha_B,... 
          alpha_B,...
          pi/3+alpha_B,... 
          pi/3-alpha_B,... 
          pi/3-alpha_B,... 
          pi/3+alpha_B];
    % gamma_P represents the direction of the points where the rod is
% attached to the platform.
    gamma_P= [pi/3-alpha_P,... 
          pi/3-alpha_P,...
          pi/3+alpha_P,... 
          alpha_P,... 
          alpha_P,... 
          pi/3+alpha_P];
      
      % the base attachment points base_i projection of joint center on corresponding axis of rotation
%       servo_attach_B= r_B* [ [cos(gamma_B(1)), -sin(gamma_B(1)), 0]',...
%                      [cos(gamma_B(2)), sin(gamma_B(2)), 0]',...
%                      [-cos(gamma_B(3)), sin(gamma_B(3)), 0]',...
%                      [-cos(gamma_B(4)), sin(gamma_B(4)), 0]',...
%                      [-cos(gamma_B(5)), -sin(gamma_B(5)), 0]',...
%                      [-cos(gamma_B(6)), -sin(gamma_B(6)), 0]'];
     servo_attach_B= r_base* [ [cos(gamma_B(1)), -sin(gamma_B(1)), 0]',...
                     [cos(gamma_B(2)), sin(gamma_B(2)), 0]',...
                     [-cos(gamma_B(3)), sin(gamma_B(3)), 0]',...
                     [-cos(gamma_B(4)), sin(gamma_B(4)), 0]',...
                     [-cos(gamma_B(5)), -sin(gamma_B(5)), 0]',...
                     [-cos(gamma_B(6)), -sin(gamma_B(6)), 0]'];
                
     % The platform attachment points pi are coincident with the centers of the corresponding joints attached to the platform       
     rod_attach_P= r_platform*[ [cos(gamma_P(1)), -sin(gamma_P(1)), 0]',...
                      [cos(gamma_P(2)), sin(gamma_P(2)), 0]',...
                      [cos(gamma_P(3)), sin(gamma_P(3)), 0]',...
                      [-cos(gamma_P(4)), sin(gamma_P(4)), 0]',...
                      [-cos(gamma_P(5)), -sin(gamma_P(5)), 0]',...
                      [cos(gamma_P(6)), -sin(gamma_P(6)), 0]'];
                  
     h= sqrt(rod_l.^2+ servo_arm_l.^2 -(rod_attach_P(1,:) - servo_attach_B(1,:)).^2 - (rod_attach_P(2,:) - servo_attach_B(2,:)).^2) -rod_attach_P(3,:);
home_pos= [0, 0, h(1)]';


%% Calculate the needed leg length
% Calculate the transformation matrix for a transform from platform to base
% system with the given euler angles. If errors occure, you have to define
% the transformation matrices around x-, y- and z-axis for yourself.
T_BP= rotZ(orient(3))*rotY(orient(2))*rotX(orient(1)); 


% Calculate the leg vector and leg length for the new position of the
% platform for each servo.
for i=1:6
     
    leg(:,i)= trans + home_pos + T_BP*rod_attach_P(:,i) - servo_attach_B(:,i);
    
    leg_length(i)= norm(leg(:,i));
end

%% Calculate the new servo angles
% Get coordinates of the points where the rod is attached to the platform 
x_P= leg(1,:) + servo_attach_B(1,:);
y_P= leg(2,:) + servo_attach_B(2,:);
z_P= leg(3,:) + servo_attach_B(3,:);

% Get coordinates of the points where the servo arm is attached to the
% servo axis.
x_B= servo_attach_B(1,:);
y_B= servo_attach_B(2,:);
z_B= servo_attach_B(3,:);

% Calculate auxiliary quatities L, N and M
L= leg_length.^2 -rod_l.^2 + servo_arm_l.^2;

M= 2*servo_arm_l*(z_P  - z_B);

for i= 1:6
    N= 2*servo_arm_l*(cos(beta(i)).*(x_P(i) - x_B(i)) + sin(beta(i)).*(y_P(i) - y_B(i)));
    
    % The wanted position could be achieved if the solution of this
    % equation is real for all i
    
    servo_angle_del(i)= asin(L(i)./sqrt(M(i)^2 + N^2))- atan2(N,M(i));
    
    % Get postion of the point where a spherical joint connects servo arm and
    % rod.

    joint_B(:,i)= [ servo_arm_l*cos(servo_angle_del(i)).*cos(beta(i)) + servo_attach_B(1,i);...
                    servo_arm_l*cos(servo_angle_del(i)).*sin(beta(i)) + servo_attach_B(2,i);...
                    servo_arm_l*sin(servo_angle_del(i))];
end
%% Plot the stewart platform

cla;
% shading interp
cmap = zeros(10, 3);
cmap = [0, 0, 1; ...   % Blue for 1
  0, 0, 0; ...       % Black for 2
  0, 1, 1; ...       % Green for 3
  0.5, 0.5, 0.5; ... % Gray for 4
  0.5, 0.5, 0; ...   % Brown for 5
  0, 0, 0; ...       % Black for 6
  0, 0, 0; ...       % Black for 7
  0, 0, 0; ...       % Black for 8
  0, 0, 0; ...       % Black for 9
  1, 0, 0]     ;      % Red for 10
colormap(cmap(10,:));

    
leg= leg + servo_attach_B;
 fill3(servo_attach_B(1,:),servo_attach_B(2,:),servo_attach_B(3,:),'-k');
 shadowcolor=.95*[1 1 1];
% patch(servo_attach_B(1,:),servo_attach_B(2,:),shadowcolor,'edgecolor','none')
 alpha(.5)
hold on
grid on
% fill3(leg(1,:),leg(2,:),leg(3,:),'-g');
% material shiny
%Lighting
%     light('style','local','position',[0 0 3],'color',[1 1 1])
    lighting gouraud
fill3(leg(1,:),leg(2,:),leg(3,:),'-r');
%  alpha(.5)
axis([ -r_base-servo_arm_l, r_base+servo_arm_l,...
       -r_base-servo_arm_l, r_base+servo_arm_l,...
       -servo_arm_l rod_l+servo_arm_l]);
rotate3d on;



for i=1:6
    line([servo_attach_B(1,i) joint_B(1,i)],... 
         [servo_attach_B(2,i) joint_B(2,i)],...
         [servo_attach_B(3,i) joint_B(3,i)],...
         'Color','r','LineWidth',3);
    plot3(joint_B(1,i),joint_B(2,i),joint_B(3,i),':ok','linewidth',2)
%     line([joint_B(1,i) leg(1,i)],... 
%          [joint_B(2,i) leg(2,i)],...
%          [joint_B(3,i) leg(3,i)],...
%          'Color','k','LineWidth',4);
     line([joint_B(1,i) leg(1,i)],... 
         [joint_B(2,i) leg(2,i)],...
         [joint_B(3,i) leg(3,i)],...
         'Color','r','LineWidth',3);
     plot3(leg(1,i),leg(2,i),leg(3,i),':ok','linewidth',2)
%     line([servo_attach_B(1,i) leg(1,i)],... 
%          [servo_attach_B(2,i) leg(2,i)],...
%          [servo_attach_B(3,i) leg(3,i)],...
%          'Color','y');
end

end