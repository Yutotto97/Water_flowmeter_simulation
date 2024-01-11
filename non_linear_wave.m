%% Define k-space
% f_max_x = medium.sound_speed/(2*dx)
% So for f = 3MHz, dx_min = 2.467e-4

% create the computational grid
Nx = 20000;       % number of grid points in the x (row) direction
dx = 0.005e-3;   % grid point spacing in the x direction [m]
kgrid = kWaveGrid(Nx, dx);

% define the properties of the propagation medium
medium.sound_speed = 1480;    % [m/s]
medium.density = 1000;        % [kg/m^3]

% Absorption and non linearity : 
% Ref : Modelling Nonlinear Ultrasound Propagation in Absorbing Media using the k-Wave Toolbox: Experimental Validation
medium.BonA = 4.96;            %non linearity parameter
medium.alpha_coeff = 2.17e-3;  % [dB/(MHz^y cm)]
medium.alpha_power = 2;   % y 

%% 2 cycles in water (freq = 3MHz)

% create initial pressure distribution using a smoothly shaped sinusoid
fsign = 6e6; %input frequency
lambda = medium.sound_speed(1)/fsign; %wavelength
width = round(lambda/dx); % [grid points]  
height = 1e6;     % [Pa]
in = (0:pi/(width/4):4*pi).';
source.p0 = [ height * sin(in) ; zeros(Nx  - width - 1, 1)];

% set the simulation time
t_end = kgrid.x_size / medium.sound_speed;

%Define time only for second order systems response
time_in = linspace(0,t_end,length(source.p0)).';

xie = 0.5;   % amortissement
w0e = 6e6;   % frequence de resonance
numerator = 1;
denominator = [1/w0e^2,2*xie/w0e,1];
syse = tf(numerator,denominator);

source.p0 = lsim(syse, source.p0, time_in); %signal acoustique a l'emission

figure;
plot(time_in, source.p0)

% create a Cartesian sensor mask
sensor.mask = [-47e-3,0, 13.6e-3];  % [mm]

% define the time array
kgrid.makeTime(medium.sound_speed*ones(Nx,1), [], t_end);

input_args = {'RecordMovie', false,'PlotPML', false 'MovieName', 'Linearity_3MHz_alpha=0.1_y=1.1_Pressure=0.1MPa', 'MovieProfile', 'MPEG-4', 'MovieArgs', {'FrameRate', 10}};
% run the simulation
sensor_data = kspaceFirstOrder1D(kgrid, medium, source, sensor, input_args{:});

figure;
plot(kgrid.t_array,sensor_data(1,:))
hold on
plot(kgrid.t_array,sensor_data(2,:))
hold on
plot(kgrid.t_array,sensor_data(3,:))
legend('NF','D=50mm','FF')
xlabel('time');
ylabel('pressure(MPa)')
%savefig('Linearity_3MHz_y=1.1.fig')

t_1=(Nx*dx/2+sensor.mask(1))/medium.sound_speed;
t_2=(Nx*dx/2+sensor.mask(2))/medium.sound_speed;
t_3=(Nx*dx/2+sensor.mask(3))/medium.sound_speed;
indice_t_1=round(t_1/kgrid.dt);
indice_t_2=round(t_2/kgrid.dt);
indice_t_3=round(t_3/kgrid.dt);

% figure;
% concatenated_output=[sensor_data(1,indice_t_1-250:indice_t_1+200), ...
%     sensor_data(2,indice_t_2-250:indice_t_2+200),sensor_data(3,indice_t_3-250:indice_t_3+200)];
% plot(concatenated_output)

figure;
resize_left=0;
resize_right=1100;
plot(kgrid.t_array(1:resize_left+resize_right+1) ...
    ,flip(sensor_data(1,indice_t_1-resize_right:indice_t_1+resize_left))*1e-6);
hold on
% plot(kgrid.t_array(1:5301),flip(sensor_data(2,indice_t_2-5000:indice_t_2+300)));
% hold on
plot(kgrid.t_array(1:resize_left+resize_right+1) ...
    ,flip(sensor_data(3,indice_t_3-resize_right:indice_t_3+resize_left))*1e-6);
legend('NF','FF');
xlabel('time');
ylabel('pressure(MPa)');
title('Hydrophone signals NF and FF (simulated)')
%savefig('Linearity_3MHz_cMUT_sensors_flipped.fig')

%% Signal input in water

Ncycle = 2;
signal_freq = 3e6;
fs = signal_freq*100;

input_signal = toneBurst(fs, signal_freq,Ncycle,'Envelope','Rectangular');

source_strength = 0.3e6;          % [Pa]

% scale the source magnitude by the source_strength divided by the
% impedance (the source is assigned to the particle velocity)
input_signal = (source_strength ./ (medium.sound_speed * medium.density)) .* input_signal;


%% Defining Transducer

Nx = 100;      % number of grid points in the x (row) direction
Ny = 1250;
Nz = 1250;
dx = 10e-6;   % grid point spacing in the x direction [m]
dy = 10e-6;
dz = 10e-6;
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

% physical properties of the transducer
transducer.number_elements = 16;    % total number of transducer elements
transducer.element_width = 374;       % width of each element [grid points]
transducer.element_length = 1200;     % length of each element [grid points]
transducer.element_spacing = 1;     % spacing (kerf width) between the elements [grid points]
transducer.radius = inf;            % radius of curvature of the transducer [m]

% calculate the width of the transducer in grid points
transducer_width = transducer.number_elements * transducer.element_width ...
    + (transducer.number_elements - 1) * transducer.element_spacing;

% use this to position the transducer in the middle of the computational grid
transducer.position = round([1, Ny/2 - transducer_width/2, Nz/2 - transducer.element_length/2]);

% properties used to derive the beamforming delays
transducer.sound_speed = 1480;                  % sound speed [m/s]
%transducer.focus_distance = 20e-3;              % focus distance [m]
%transducer.elevation_focus_distance = 19e-3;    % focus distance in the elevation plane [m]
%transducer.steering_angle = 0;                  % steering angle [degrees]

% % define the transducer elements that are currently active
% transducer.active_elements = zeros(transducer.number_elements, 1);
% transducer.active_elements(:) = 1;

% append input signal used to drive the transducer
transducer.input_signal = input_signal;

% Define time array
t_end = kgrid.x_size / medium.sound_speed;  
dt = 2e-9;  % [s]
kgrid.setTime(round(t_end / dt) + 1, dt);

% create the transducer using the defined settings
transducer = kWaveTransducer(kgrid, transducer);

y = (-round(Ny/4):2:round(Ny/4)) * dy;          % [m]
z = (Nz-1)*dz/2*ones(length(y));  % [m]
x = (Nx-1)*dx/2*ones(length(z));  % [m]
sensor.mask = [x; y; z];

%[grid_data, order_index, reorder_index] = cart2grid(kgrid, cart_data)

input_args = {'PlotLayout', true,'DataCast', 'single','PlotPML',false,'RecordMovie', true, 'MovieName', '3D_transducer_example', 'MovieProfile', 'MPEG-4', 'MovieArgs', {'FrameRate', 10}};
% run the simulation
sensor_data = kspaceFirstOrder3D(kgrid, medium, transducer,sensor,input_args{:});

