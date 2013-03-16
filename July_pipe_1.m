% Function to calculate volumetric surface strain from a
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

clear all
close all
%Inputs

first_y_values = [];  %empty value to establish existence
second_y_values = [];  %empty value to establish existence
second_nu_values = [];
second_p_values = [];
first_nu_values = [];
first_p_values = [];

hold on
	
c1 = 10;
c2 = 1000;
x = -6500:5:6500;
y = -6500:5:6500;
%r = -50000:50:50000;
r = sqrt(x.^2 + y.^2);
R1 = sqrt(r.^2 + c1.^2);
R2 = sqrt(r.^2 + c2.^2); 
%a = 15;
P = 600e6:10e6:910e6;
nu = 0.24:0.005:0.33;
G = 5e9;
experiment_number = 1;

%---------------------------------------------------------------------------------------------------------
%beginning of loop for all values of nu_index and p_index

for p_index=1:max(size(P)) %for all values of P from the first value to the maximum (last value)
  for nu_index=1:max(size(nu)) %for all values of nu from the first value to the maximum (last value)
    figure(experiment_number)  % show the resultant graph for each experiment sequentially
	hold on
	
	%----------------------------------------------------
	%calculations for ur
    a = 15;							%running a loop, need to respecify within loop or a will become a matrix
	ur = (a^2*P(p_index))/(4*G);
	j = (c1.^3)./R1.^3; 
	k = (2*c1*(-3+5*nu(nu_index)))./R1; 
	l = (5*c2.^3*(1-2*nu(nu_index))-2*c2*r.^2*(-3+5*nu(nu_index)))./R2.^3; 
	ur = ur*(j+k+l);
	ur = ur.*(x./r.^2);


	%------------------------------------------------------
	%calculations for volumetric pipe strain

	 
	%E00 = ur/x
	tang_strain=ur./x;

	%Err = diff ur/ diff x

	a = diff(ur);
	b = diff(x);

	rad_strain=a./b;

	%EV = (1-2*nu/1-nu)*(Err+E00)
	volumetric_pipe=(1-2*nu(nu_index)./1-nu(nu_index))*(rad_strain + tang_strain(:,1:max(size(tang_strain))-1));

	vent_distance = x;    %distance from the vent 

	plot (vent_distance(:,1:max(size(vent_distance))-1),volumetric_pipe)
	title_label = sprintf('p=%f nu=%f', P(p_index), nu(nu_index));    %inserting title on graphs for values of p and nu used
	title(title_label)

	for i=1:max(size(vent_distance))
	  if (vent_distance(i) == 4600)            %for x = 4600, use the accompanying y value in the operations below
	    fprintf('4600: nu=%f, p=%d, y=%d\n', nu(nu_index), P(p_index), volumetric_pipe(i)); %print the following values as the loop is run
		first_y_values = [first_y_values volumetric_pipe(i)];
		first_nu_values = [first_nu_values nu(nu_index)];
		first_p_values = [first_p_values P(p_index)];
	    if ((volumetric_pipe(i) > 4e-7) && (volumetric_pipe(i) < 4e-9))
			scatter(vent_distance(i), volumetric_pipe(i))  %if the y value produced is in the range specified plot its coordinates
			lbl = sprintf('x=%d, y=%d, nu=%f, p=%f', vent_distance(i), volumetric_pipe(i), nu(nu_index), P(p_index))  
			text(vent_distance(i)-4000, volumetric_pipe(i)+(0.2e-4), lbl) % text label placed at these positions
			vent_distance(i)
			volumetric_pipe(i)
		end % if y is in range
	  end
	  if (vent_distance(i) == 6000)
	    fprintf('6000: nu=%f, p=%d, y=%d\n', nu(nu_index), P(p_index), volumetric_pipe(i));
		second_y_values = [second_y_values volumetric_pipe(i)];
		second_nu_values = [second_nu_values nu(nu_index)];
		second_p_values = [second_p_values P(p_index)];
			if ((volumetric_pipe(i) > 4e-7) && (volumetric_pipe(i) < 4e-9))
			scatter(vent_distance(i), volumetric_pipe(i))
			lbl = sprintf('x=%d, y=%d, nu=%f, p=%f', vent_distance(i), volumetric_pipe(i), nu(nu_index), P(p_index))
			text(vent_distance(i)-4000, volumetric_pipe(i)-(0.2e-4), lbl)
			vent_distance(i)
			volumetric_pipe(i)
		end % if y is in range
	  end
	end 
	close(experiment_number) % close each experiment graph before the next one opens 
	experiment_number = experiment_number + 1;
  end % nu iterator end
end % p iterator end

% Get experiment ID of the minimum value at 4600 and print
min_exp_id = find(first_y_values == min(first_y_values));
fprintf('Min at 4600 is where nu=%f, P=%d and y=%d\n', first_nu_values(min_exp_id), first_p_values(min_exp_id), first_y_values(min_exp_id))

% Get experiment ID of the minimum value at 6000 and print
max_exp_id = find(second_y_values == min(second_y_values));
fprintf('Min at 6000 is where nu=%f, P=%d and y=%d\n', second_nu_values(max_exp_id), second_p_values(max_exp_id), second_y_values(max_exp_id))

fprintf('P was in the range %d to %d\n', min(P), max(P))
fprintf('Nu was in the range %f to %f\n', min(nu), max(nu))
% Print information about minimum 4600
%fprintf('Min at 4600 is %d\n', min(first_y_values))   % show min y value for all calculated at 4600m
%fprintf('Min at 6000 is %d\n', min(second_y_values))   % show min y value for all calculated at 6000m

	fprintf('4600) chi squared = %d\n', (-3.5e-8 - first_y_values(min_exp_id))^2/-first_y_values(min_exp_id));
	fprintf('6000) chi squared = %d\n', (-1.2e-8 - first_y_values(min_exp_id))^2/-first_y_values(min_exp_id));
	chi_4600=(-3.5e-8 - first_y_values(min_exp_id))^2/-first_y_values(min_exp_id);
	chi_6000= (-1.4e-8 - first_y_values(min_exp_id))^2/-first_y_values(min_exp_id);
	chi_4600+chi_6000
	