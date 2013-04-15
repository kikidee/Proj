% Function to calculate the tangential surface strain 
% from a spherical source at depth

%--------------------------------------------------

%variables:

% d = source depth [m]
% a = source radius [m]
% P = source pressure difference from surround [Pa]
% nu = Poisson Ratio
% G = Shear Modulus (Rigidity) [Pa]
% E = Young's Modulus [Pa]
% x = radial distance on the surface [m]
% ur = radial displacement

%----------------------------------------------------

clear all
close all
%Inputs

d = 5000;
a = 1000;    
P = 10e6;
nu = 0.25;
E = 20e9;
G = 8e9;
x = -15000:100:15000;
y = -15000:100:15000;

c = x;
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

%calcuate ur from mogi equation
ur = P.*(a.^3).*(1-nu)./G.*(c./(c.^2+d.^2).^1.5);



% Calculations

% E00 = radial displacement/radial distance
Tstrain = ur./c;

% normalising factor
Normal = (3*P*(a^3))/(4*G*d^3);

norm_tang_strain = Tstrain/((3*P*a^3)/(4*G*d^3));

source_depth = c/d;

plot(source_depth,norm_tang_strain, 'LineWidth', 2)

%----------------------------------------------------------------

title('Tangential Strain (Lisowski)', 'FontSize', 12, 'FontName', 'Arial');		% Set Graph Title in fontsize

ylabel('Tangential Strain (strain units)', 'FontSize', 12)						% Set Y Axis label
xlabel('Distance (Top Depths)', 'FontSize', 12)										% Set X Axis label

% Set Graph Limits
%xMin xMax yMin yMax
%axis([3 3 -1 12])


waitforbuttonpress()

%Save figure 1 to jpeg,called lisowkski_tang_strain.jpg
% at a resolution of 500 dots per inch
% text is (for commercial print) 300
% images are 2000
print(1,'-djpeg','lisowkski_tang_strain','-r500')





