%% simulate in 2D

% create the computational grid
Nx = 1500;           % number of grid points in the x (row) direction
Ny = 1500;           % number of grid points in the y (column) direction
dx = 10e-6;        % grid point spacing in the x direction [m]
dy = 10e-6;        % grid point spacing in the y direction [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy);

medium.sound_speed = 1480;    % [m/s]
medium.density = 1000;        % [kg/m^3]

% Absorption and non linearity : 
% Ref : Modelling Nonlinear Ultrasound Propagation in Absorbing Media using the k-Wave Toolbox: Experimental Validation
medium.BonA = 4.96;            %non linearity parameter
medium.alpha_coeff = 2.17e-3;  % [dB/(MHz^y cm)]
medium.alpha_power = 2;   % y 

% create initial pressure distribution using a smoothly shaped sinusoid
fsign = 3e6; %input frequency
lambda = medium.sound_speed(1)/fsign; %wavelength
width = round(lambda/dx); % [grid points]  
height = 100;     % [Pa]
in = (0:pi/(width/4):4*pi).';

R_t=Ny*3;   %résolution temporelle

source_elec = [ height * sin(in) ; zeros(R_t  - width - 1, 1)];

% set the simulation time
t_end = kgrid.y_size / medium.sound_speed;

%Define time only for second order systems response
kgrid.t_array = linspace(0,t_end,R_t).';

xie = 0.5;   % amortissement
w0e = 3e6;   % frequence de resonance
numerator = 1;
denominator = [1/w0e^2,2*xie/w0e,1];
syse = tf(numerator,denominator);

% define source mask for a linear transducer with an odd number of elements   
num_elements = 320;      % [grid points]
element_width = 3;
element_spacing = 1;
source.p_mask = zeros(Nx, Ny);
element = [ones(element_width,1);zeros(element_spacing,1)];
array = [];
for i=1:num_elements
    array =  [array;element];
end
start_index = Nx/2 - round(length(array)/2) ;
source.p_mask(start_index : start_index + length(array) -1, 1) = array;

source.p = lsim(syse, source_elec, kgrid.t_array).'; %signal acoustique a l'emission

x=zeros(1,round(Nx/100));
y=linspace(-(Ny-10)*dy/2,(Ny-10)*dy/2,length(x));
sensor.mask = [x; y];

% % each sensor point
% sensor.record = {'p_final', 'p_max', 'p_rms'};

% create a display mask to display the transducer
display_mask = source.p_mask;
 
% assign the input options
input_args = {'PMLInside', false, 'PlotPML', false,'RecordMovie', true, 'MovieName', 'cMUT_2D_wave_propagation_pa=100_d=1.5cm_side', 'MovieProfile', 'MPEG-4', 'MovieArgs', {'FrameRate', 1}};

% run the simulation
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

figure;
plot(kgrid.t_array,sensor_data(2,:))
hold on
plot(kgrid.t_array,sensor_data(11,:))
%norm_data=sensor_data.p_final/max(max(sensor_data.p_final));

% figure;
% s=surf(kgrid.x_vec,kgrid.y_vec,sensor_data.p_final.');
% s.EdgeColor='none';



%% simulate in 2D

% create the computational grid
Nx = 2*1250;           % number of grid points in the x (row) direction
Ny = 2*1000;           % number of grid points in the y (column) direction
dx = 5e-6;        % grid point spacing in the x direction [m]
dy = 5e-6;        % grid point spacing in the y direction [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy);

medium.sound_speed = 1480;    % [m/s]
medium.density = 1000;        % [kg/m^3]

% Absorption and non linearity : 
% Ref : Modelling Nonlinear Ultrasound Propagation in Absorbing Media using the k-Wave Toolbox: Experimental Validation
medium.BonA = 4.96;            %non linearity parameter
medium.alpha_coeff = 2.17e-3;  % [dB/(MHz^y cm)]
medium.alpha_power = 2;   % y 

% create initial pressure distribution using a smoothly shaped sinusoid
fsign = 3e6; %input frequency
lambda = medium.sound_speed(1)/fsign; %wavelength
width = round(lambda/dx); % [grid points]  
height = 0;     % [Pa]
in = (0:pi/(width/4):4*pi).';

R_t=Ny*3;   %résolution temporelle

source_elec = [ height * sin(in) ; zeros(R_t  - width - 1, 1)];

% set the simulation time
t_end = kgrid.y_size / medium.sound_speed;

%Define time only for second order systems response
kgrid.t_array = linspace(0,t_end,R_t).';

xie = 0.5;   % amortissement
w0e = 11e6;   % frequence de resonance
numerator = 1;
denominator = [1/w0e^2,2*xie/w0e,1];
syse = tf(numerator,denominator);

% define source mask for a linear transducer with an odd number of elements   
num_elements = 16;      % [grid points]
element_width = 20;
element_spacing = 1;
source.p_mask = zeros(Nx, Ny);
element = [ones(element_width,1);zeros(element_spacing,1)];
array = [];
for i=1:num_elements
    array =  [array;element];
end
start_index = Nx/2 - round(length(array)/2) ;
source.p_mask(start_index : start_index + length(array) -1, 1) = array;

source.p = lsim(syse, source_elec, kgrid.t_array).'; %signal acoustique a l'emission

x=zeros(1,round(Nx/100));
y=linspace(-(Ny-10)*dy/2,(Ny-10)*dy/2,length(x));
sensor.mask = [x; y];

% each sensor point
sensor.record = {'p_final', 'p_max', 'p_rms'};

% create a display mask to display the transducer
display_mask = source.p_mask;

% assign the input options
input_args = {'DisplayMask', display_mask, 'PMLInside', false, 'PlotPML', false,'RecordMovie', true, 'MovieName', 'cMUT_2D_wave_propagation_pa=100000', 'MovieProfile', 'MPEG-4', 'MovieArgs', {'FrameRate', 10}};

% run the simulation
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

% figure;
% plot(kgrid.t_array,sensor_data(2,:))

%norm_data=sensor_data.p_final/max(max(sensor_data.p_final));

% figure;
% surf(kgrid.x_vec,kgrid.y_vec,sensor_data.p_final.')
% colorbar
