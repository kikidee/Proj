% Function to calculate tilt from a
% prolate source at depth

%--------------------------------------------------

%variables:

%c1 = top depth of pipe
%c2 = bottom of the pipe
%r = sqrt(x^2 + y^2)
%R1 = sqrt(r^2 + c1^2)
%R2 = sqrt(r^2 + c2^2) 
% a = source radius [m]
% P = source pressure difference from surround [Pa]
% nu = Poisson Ratio
% G = Shear Modulus (Rigidity) [Pa]
% E = Young's Modulus [Pa]
% x = radial distance on the surface [m]
% y = ?
% ur = radial displacement

%----------------------------------------------------

% Tidy Up
%--------
clear all
close all

% Constants
c1 = 1000;
c2 = 50000;
a = 500;
P = 10e6;
nu = [0.25 0.225 0.275];
G = 8e9;

% X & Y Displacement Vectors
x = -4000:5:4000;
y = -4000:5:4000;
close all
figure(1)

hold on
for nu_value=1:max(size(nu))
	% Distance Vector (+ AND -)
	for i=1:max(size(x))
	  if (x(i) < 0) % dont care about y as y=x
		r(i)=-sqrt(x(i)^2 + y(i)^2);
	  else
		r(i)=sqrt(x(i)^2 + y(i)^2);
	  end
	end

	R1 = sqrt(r.^2 + c1.^2);
	R2 = sqrt(r.^2 + c2.^2); 

	%----------------------------------------------------
	%calculations for w
	w = -(a^2*P)/(4*G);
	j = (c1^2)./R1.^3; 
	k = (2*(-2+5*nu(nu_value)))./R1; 
	l = (c2^2*(3-10*nu(nu_value))-2*r.^2*(-2+5*nu(nu_value)))/R2.^3; 
	w = w*(j+k+l);

	%-------------------------------------------------------
	%calculations for tilt (wx,wy)

	top_depth = -r/c1;
	  
	alpha = diff(w);
	beta = diff(x);
	gamma = diff(y);
	tilt_x = -(alpha./beta);
	tilt_y = -(alpha./gamma);

	% Remove 1 value of top_depth due to diff function
	top_depth = top_depth(:,1:max(size(top_depth))-1);

	axis([-4 4 -5 5])
	%aacolor=('blue', 'blue', 'blue');
	%printf('color is %d\n', aacolor(nu_value));
	if (nu_value == 1)
	  plot(top_depth, tilt_x*100000, 'blue', 'LineWidth',2)
	elseif (nu_value == 2)
	  plot(top_depth, tilt_x*100000, 'red', 'LineWidth',2)
	else
	  plot(top_depth, tilt_x*100000, 'green', 'LineWidth',2)
	end
end

% Set Graph Title in fontsize
title('Tilt (Lisowski)', 'FontSize', 12, 'FontName', 'Arial');

legend('Tilt (nu = 0.25)','Tilt (nu = 0.225)','Tilt (nu = 0.275)')
xlabel('Distance (Top Depths)', 'FontSize', 12, 'FontName', 'Arial')
ylabel('Tilt', 'FontSize', 12, 'FontName', 'Arial')

waitforbuttonpress()

%----------------------------------------------------------------------------------
% Save high resolution version of graph to working directory
%----------------------------------------------------------------------------------

% Save figure 1 to jpeg,called output_test.jpg
% at a resolution of 500 dots per inch
% text is (for commercial print) 300
% images are 2000
print(1,'-djpeg','lisowski_pipe_tilt','-r500')

%plot(top_depth(:,1:max(size(top_depth))-1), tilt_y_a)

%grtoup.function(argument)
%a= diff.(w);
%b= diff.(x);

%tilt_x = -(diff.(w)/diff.(x));
