%===================================================================================================================
% Adapted script from the Yang Model, written to produce 2D graphs of displacement and strain produced by 
% inflation/deflation of a shallowly dipping finite prolate spheroid and directly compare them with COMSOL output
%===================================================================================================================

clear all;
close all;

%===================================================================================================================
% Inputs
%===================================================================================================================

eps=1e-10;   			% small number
z0=.51;                	% center depth (positive)
P= 12900;               %a value of 500 gives 25MPa for mu=5GPa, for 10 GPa divide value by 10 and result is MPa for 5GPa divide by 20
							% excess pressure, mu*10^(-5) Pa
a =0.5;                	% major axis, km
b=.015;                 % minor axis, km
theta=pi/180*89.999;  	% error in solution for 90 degrees use 89.9 instead for vertical prolate source
						% strike, rad  (0-2*pi)
                      
%===================================================================================================================
% Axis Settings
%===================================================================================================================

x_step = .005; % Resolution of X Dimension in KM
y_step = .005; % Resolution of Y Dimension in KM
x_point_count = 201; % Number of points in X Dimension
y_point_count = 201; % Number of points in Y Dimenion

%===================================================================================================================
% Generate 3D Grid From Axis
%===================================================================================================================

% Build From 0 -> N then shift by -N/2
N = (y_point_count-1)*x_step;

%                -N/2         step      +N/2
x_data_points = -(N/2) : x_step : (N/2);
[x,y] = meshgrid(x_data_points,x_data_points');
flat_z_dimension = zeros(size(x)); % Generate All-Zero Matrix of Z-Values

gamma = 1;    % 1st Lame constant
mu    = 1;    % shear modulus (2nd Lame constant)
nu    = 0.25; % Poisson's ratio

coeffs(1)=1/(16*mu*(1-nu));
coeffs(2)=3-4*nu;
coeffs(3)=4*(1-nu)*(1-2*nu);

c=sqrt(a^2-b^2);
[sph]=spheroid(a,b,c,[gamma mu nu],0,theta,P); % 0 was Phi

[Up1,Up2,Up3] = yang(sph,c,z0,x,x',0,[gamma mu nu],[sin(theta) cos(theta)],coeffs,flat_z_dimension);
[Um1,Um2,Um3] = yang(sph,-c,z0,x,x',0,[gamma mu nu],[sin(theta) cos(theta)],coeffs,flat_z_dimension);

ux = -double(Up1) + double(Um1);  	%displacement in the x direction
uz =  Up3 - Um3;					%displacement in the z direction


%tang_yang = tang_disp(101,:)./x(101,:);
%rad_yang=diff(tang_disp(101,101:201))./diff(x(101,101:201));

%vol_yang = (1-2*nu/1-nu)*(tang_yang(101:200)+rad_yang);

%=====================================================================================================================
% Calculate Surface Displacement
%=====================================================================================================================

tang_disp = sqrt(ux.^2+(ux').^2); 									%inserted by jo to quantify horizontal displacement
%vertyang       = uz(101,:)/x(101,:);                              	% vertical displacement from Yang Model u
%vert_yang      = (nu/(1-nu))*(rad_yang+tang_yang(101:200));			% vertical displacement from strain components (Mogi Model)

%=====================================================================================================================
% Call Comsol data (positive pressure)
%=====================================================================================================================

alldatarray    = load('Pipe_correct_depth_by_hand.txt','-ascii');             	% exported data for positive pressure (645.5Mpa)
dist           = alldatarray(:,1);
ustrain        = alldatarray(:,2);
wstrain        = alldatarray(:,3);
volstrain      = alldatarray(:,4);

%=====================================================================================================================
% Call Comsol data (negative pressure)
%=====================================================================================================================

%alldatarrayneg = load('july_bench_negative.txt','-ascii');   		% exported data for negative pressure (-645.5Mpa)
%distneg        = alldatarrayneg(:,1);
%ustrainneg     = alldatarrayneg(:,3);
%wstrain       = alldatarray(:,4);
%volstrainneg   = alldatarrayneg(:,2);

%=====================================================================================================================
% Establish Strain Equations
%=====================================================================================================================

tang_yang      = tang_disp(101,:)./x(101,:);						% tangential strain equation from Mogi Model
rad_yang       = diff(tang_disp(101,101:201))./diff(x(101,101:201));% radial strain equation from Mogi Model
vol_yang       = (1-2*nu/1-nu)*(tang_yang(101:200)+rad_yang);		% volumetric strain equation from Bonaccorso Model
%vol_ymogi      = (tang_yang(101:200)+rad_yang+vert_yang);          % mogi equation for volumetric strain
%vol_yangedit   = (tang_yang(101:200)+rad_yang);					% edited yang equation 


tang_com       = ustrain./(dist);									% tangential strain equation from Mogi Model
rad_com        = diff(ustrain)./diff((dist));						% radial strain equation from Mogi Model
vol_com        = ((1-2*nu)/(1-nu))*(tang_com(1:200)+rad_com);			% volumetric strain equation from Bonaccorso Model

%tang_comneg    = ustrainneg*100./(distneg/1000);
%rad_comneg     = diff(ustrainneg*100)./diff((distneg/1000));
%vol_comneg     = (1-2*nu/1-nu)*(tang_comneg(1:100)+rad_comneg);

%===================================================================================================================
% Plot Displacement and Strain
%===================================================================================================================

plot(x(101,101:200)*1000, vol_yang/100, 'LineWidth',2)
hold on
%x_fig = 0;
%y_fig = 0;
%spacing = 0.001;
%for i=1:max(size(rad_com))											% scatter comsol data evenly spaced
%  x_new = dist(i)/max(dist);										
%  y_new = rad_com(i)/(2*max(rad_com));								% Y axis data is very small, therefore its weight needs to be increased
%  if (abs(x_new-x_fig)^2 + abs(y_new-y_fig)^2) > spacing
%    x_fig = x_new;
%	y_fig = y_new;
%   scatter(dist(i),rad_com(i),50,'g', 'fill')
%  end
%end
%===================================================================================================================
% Plot Volumetric Strain
%===================================================================================================================


%plot(dist,volstrain, 'red')
plot(dist,volstrain, 'green')
%hold on
%plot(dist(1:200),vol_com, 'red')
%plot(x(101,101:200), vol_yang, 'blue')
%plot(x(101,101:200)*1000, vol_ymogi/100000, 'red')
%plot(x(101,101:200)*1000, vol_yangedit/100000, 'black')

%===================================================================================================================
% Tidy Plot
%===================================================================================================================

%title('Radial Strain at -500m Top Depth', 'FontSize', 12, 'FontName', 'Arial');	% Set Graph Title in fontsize

%ylabel('Strain (strain units)', 'FontSize', 12)									% Set Y Axis label
%xlabel('Distance (metres)', 'FontSize', 12)										% Set X Axis label

%waitforbuttonpress()

line_1_name = 'Analytical';														%legend commands
line_2_name = 'Numerical';
legend(line_1_name,line_2_name, 'Location','NorthEast')

% Save figure 1 to jpeg,called output_pipe_10_deep.jpg
% at a resolution of 500 dots per inch
% text is (for commercial print) 300
% images are 2000
%print(1,'-djpeg','output_pipe_500_deep_rad','-r500')

%======================================================================================================================

