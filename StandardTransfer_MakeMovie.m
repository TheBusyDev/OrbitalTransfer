% Create a video for the Standard Transfer
% WARNING:  this script needs 'StandardTransfer.m' to be already run,
%           so that all the necessary variables are initialized.

%% Define Transfer Data
% Orbital elements
a =     [a_i,   a_i,    a_i,    a_t3,   a_f];
e =     [e_i,   e_i,    e_i,    e_t3,   e_f];
i =     [i_i,   i_f,    i_f,    i_f,    i_f];
OM =    [OM_i,  OM_f,   OM_f,   OM_f,   OM_f];
om =    [om_i,  om_2,   om_3,   om_t3,  om_f];

% True anomaly of each maneuver
theta_man = [
    theta_i,    theta_2
    theta_2,    theta_3
    theta_4,    2*pi
    0,          pi
    pi,         theta_f
];

% color of each orbit
color = [color_i, color_2, color_3, color_4, color_f];

% True anomaly of each point to plot
theta_point(1).theta(2).f = theta_i;
theta_point(1).label(2).f = "initial point";
theta_point(1).color(2).f = color(1);

theta_point(2).theta(1).f = [theta_2a, theta_2b];
theta_point(2).label(1).f = ["point A", "point B"];
theta_point(2).color(1).f = [color_w, color_c];

theta_point(2).theta(2).f = theta_2;
theta_point(2).label(2).f = "theta 2";
theta_point(2).color(2).f = color(2);

theta_point(3).theta(1).f = [theta_3c_f, theta_3d_f];
theta_point(3).label(1).f = ["point C", "point D"];
theta_point(3).color(1).f = [color_c, color_w];

theta_point(3).theta(2).f = theta_4;
theta_point(3).label(2).f = "theta 4";
theta_point(3).color(2).f = color(3);

theta_point(4).theta(1).f = [0, pi];
theta_point(4).label(1).f = ["pericenter orbit 3", "apocenter final orbit"];
theta_point(4).color(1).f = [color_w, color_c];

theta_point(5).theta(1).f = theta_f;
theta_point(5).label(1).f = "final point";
theta_point(5).color(1).f = color(5);

% label to display in legend
orbit_label = ["Initial Orbit", "Orbit 2", "Orbit 3",...
                "Bitangent Transfer Orbit", "Final Orbit"];

% Other useful variables
dth = pi/50;            % true anomaly discretisation step
dt = 70;                % time discretisation step
color_s = "#FF7300";    % satellite color
line_width = 1;
len = length(a);        % number of orbits to plot
az = 88.067;            % azimuth angle
el = 40.887;            % elevation angle
figure_pos = [0 0 750 750];         % set figure position
axis_limit = [-2 2 -2 2 -1 1]*10^4; % set axis limits


%% Generate a video of the whole transfer
clc; close all

% Video & figure initialization
video = VideoWriter("video");
open(video)
figure ('Position', figure_pos)

% Initialize true anomaly array 'theta'
theta = [];

for j = 1:len
    theta_tmp = GetTrueAnomaly (theta_man(j, 1), theta_man(j, 2), dt, a(j), e(j), mu);
    theta = [theta, theta_tmp];

    man_point(j) = length(theta);   % index of maneuver point
end

man_point = man_point(1 : end-1);

% Plot all the orbits
for j = 1:len

    % Initial Orbit
    if j == 1
        [orbit, point] = plotOrbit(a(j), e(j), i(j), OM(j), om(j), theta_man(j, 1), theta_man(j, 2), dth, mu, orbit_label(j),...
            theta(1), "initial point");

        set (point, 'Color', color(j))
        
    % Final Orbit
    elseif j == len
        [orbit, point] = plotOrbit(a(j), e(j), i(j), OM(j), om(j), theta_man(j, 1), theta_man(j, 2), dth, mu, orbit_label(j),...
            theta(end), "final point");

        set (point, 'Color', color(j))

    % Other Transfer Orbits
    else
        orbit = plotOrbit(a(j), e(j), i(j), OM(j), om(j), theta_man(j, 1), theta_man(j, 2), dth, mu, orbit_label(j));
    end

    set (orbit, 'Color', color(j), 'LineStyle', '--', 'LineWidth', line_width)
end

% plot Earth
plotEarth (Re, 0.5);

% plot properties
view(az, el)        % adjust plot3 view
axis(axis_limit)    % set axis limits
legend off

% Plot satellite and generate a video from each plot
% Satellite Initial Orbital Elements
a_s =       a(1);
e_s =       e(1);
i_s =       i(1);
OM_s =      OM(1);
om_s =      om(1);

% Satellite orbit color
color_o =   color(1);

% Initialize satellite
rr_sat_old = parorb2rv (a_s, e_s, i_s, OM_s, om_s, theta(1), mu);    % get satellite vectorial position
sat = plot3 (0, 0, 0, 'o', 'LineWidth', 1.5, 'Color', color_s, 'HandleVisibility', 'off');

for k = 1:length(theta)
    % draw satellite on its current position
    rr_sat = parorb2rv (a_s, e_s, i_s, OM_s, om_s, theta(k), mu);    % get satellite vectorial position
    set (sat, 'XData', rr_sat(1), 'YData', rr_sat(2), 'ZData', rr_sat(3))

    % draw satellite orbit
    x = [rr_sat_old(1), rr_sat(1)];
    y = [rr_sat_old(2), rr_sat(2)];
    z = [rr_sat_old(3), rr_sat(3)];
    plot3 (x, y, z, 'Color', color_o, 'LineWidth', 1.5, 'HandleVisibility', 'off');

    rr_sat_old = rr_sat;

    % get frame for the video
    frame = getframe(gcf);
    writeVideo (video, frame);

    % set satellite orbital elements for the next iteration
    index = find( k == man_point ) + 1;

    if ~isempty (index)
        a_s =       a(index);
        e_s =       e(index);
        i_s =       i(index);
        OM_s =      OM(index);
        om_s =      om(index);
        color_o =   color(index);
    end
end

close(video)


%% Generate a video of each maneuver (except for Bitangent)
dt = 50;    % time discretisation step

for j = 1:2
    clc; close all

    % Video & figure initialization
    filename = sprintf('video_%d', j);
    video = VideoWriter(filename);
    open(video)
    figure ('Position', figure_pos)

    % Plot orbits
    % plot first orbit
    [orbit, point] = plotOrbit(a(j), e(j), i(j), OM(j), om(j), 0, 2*pi, dth, mu, orbit_label(j),...
        theta_point(j).theta(2).f, theta_point(j).label(2).f);

    % set orbit color
    set (orbit, 'Color', color(j), 'LineStyle', '--', 'LineWidth', line_width);
    
    % set point color
    for n = 1:length(point)
        set (point(n), 'Color', theta_point(j).color(2).f(n))
    end

    % plot second orbit
    [orbit, point] = plotOrbit(a(j+1), e(j+1), i(j+1), OM(j+1), om(j+1), 0, 2*pi, dth, mu, orbit_label(j+1),...
        theta_point(j+1).theta(1).f, theta_point(j+1).label(1).f);
    
    % set orbit color
    set (orbit, 'Color', color(j+1), 'LineStyle', '--', 'LineWidth', line_width);

    % set point color
    for n = 1:length(point)
        set (point(n), 'Color', theta_point(j+1).color(1).f(n))
    end

    % plot properties
    view(az, el)        % adjust plot3 view
    axis(axis_limit)    % set axis limits

    % calculate true anomaly 'theta'
    theta = GetTrueAnomaly (theta_man(j, 1), theta_man(j, 2), dt, a(j), e(j), mu);

    % Satellite Initial Orbital Elements
    a_s =       a(j);
    e_s =       e(j);
    i_s =       i(j);
    OM_s =      OM(j);
    om_s =      om(j);

    % Satellite orbit color
    color_o =   color(j);

    % Initialize satellite
    rr_sat_old = parorb2rv (a_s, e_s, i_s, OM_s, om_s, theta(1), mu);    % get satellite vectorial position
    sat = plot3 (0, 0, 0, 'o', 'LineWidth', 1.5, 'Color', color_s, 'HandleVisibility', 'off');

    for k = 1:length(theta)
        % draw satellite on its current position
        rr_sat = parorb2rv (a_s, e_s, i_s, OM_s, om_s, theta(k), mu);    % get satellite vectorial position
        set (sat, 'XData', rr_sat(1), 'YData', rr_sat(2), 'ZData', rr_sat(3))

        % draw satellite orbit
        x = [rr_sat_old(1), rr_sat(1)];
        y = [rr_sat_old(2), rr_sat(2)];
        z = [rr_sat_old(3), rr_sat(3)];
        plot3 (x, y, z, 'Color', color_o, 'LineWidth', 1.5, 'HandleVisibility', 'off');

        rr_sat_old = rr_sat;

        % get frame for the video
        frame = getframe(gcf);
        writeVideo (video, frame);
    end

    close(video)
end


%% Generate a video of Bitangent maneuver
clc; close all

dt = 70;    % time discretisation step
j = 3;      % orbit index

% Video & figure initialization
filename = sprintf('video_%d', j);
video = VideoWriter(filename);
open(video)
figure ('Position', figure_pos)

% Plot orbits
% plot first orbit
[orbit, point] = plotOrbit(a(j), e(j), i(j), OM(j), om(j), 0, 2*pi, dth, mu, orbit_label(j),...
    theta_point(j).theta(2).f, theta_point(j).label(2).f);

% set orbit color
set (orbit, 'Color', color(j), 'LineStyle', '--', 'LineWidth', line_width);

% set point color
for n = 1:length(point)
    set (point(n), 'Color', theta_point(j).color(2).f(n))
end

% plot second orbit (Bitangent Transfer)
[orbit, point] = plotOrbit(a(j+1), e(j+1), i(j+1), OM(j+1), om(j+1), theta_man(j+1, 1), theta_man(j+1, 2), dth, mu, orbit_label(j+1),...
    theta_point(j+1).theta(1).f, theta_point(j+1).label(1).f);

% set orbit color
set (orbit, 'Color', color(j+1), 'LineStyle', '--', 'LineWidth', line_width);

% set point color
for n = 1:length(point)
    set (point(n), 'Color', theta_point(j+1).color(1).f(n))
end

% plot third orbit
[orbit, point] = plotOrbit(a(j+2), e(j+2), i(j+2), OM(j+2), om(j+2), 0, 2*pi, dth, mu, orbit_label(j+2),...
    theta_point(j+2).theta(1).f, theta_point(j+2).label(1).f);

% set orbit color
set (orbit, 'Color', color(j+2), 'LineStyle', '--', 'LineWidth', line_width);

% set point color
for n = 1:length(point)
    set (point(n), 'Color', theta_point(j+2).color(1).f(n))
end

% plot properties
view(az, el)        % adjust plot3 view
axis(axis_limit)    % set axis limits

% calculate true anomaly 'theta'
theta = GetTrueAnomaly (theta_man(j, 1), theta_man(j, 2), dt, a(j), e(j), mu);
man_point = length (theta);

theta_tmp = GetTrueAnomaly (theta_man(j+1, 1), theta_man(j+1, 2), dt, a(j+1), e(j+1), mu);
theta = [theta, theta_tmp];
man_point(2) = length (theta);

theta_tmp = GetTrueAnomaly (theta_man(j+2, 1), theta_man(j+2, 2), dt, a(j+2), e(j+2), mu);
theta = [theta, theta_tmp];

% Satellite Initial Orbital Elements
a_s =       a(j);
e_s =       e(j);
i_s =       i(j);
OM_s =      OM(j);
om_s =      om(j);

% Satellite orbit color
color_o =   color(j);

% Initialize satellite
rr_sat_old = parorb2rv (a_s, e_s, i_s, OM_s, om_s, theta(1), mu);    % get satellite vectorial position
sat = plot3 (0, 0, 0, 'o', 'LineWidth', 1.5, 'Color', color_s, 'HandleVisibility', 'off');

for k = 1:length(theta)
    % draw satellite on its current position
    rr_sat = parorb2rv (a_s, e_s, i_s, OM_s, om_s, theta(k), mu);    % get satellite vectorial position
    set (sat, 'XData', rr_sat(1), 'YData', rr_sat(2), 'ZData', rr_sat(3))

    % draw satellite orbit
    x = [rr_sat_old(1), rr_sat(1)];
    y = [rr_sat_old(2), rr_sat(2)];
    z = [rr_sat_old(3), rr_sat(3)];
    plot3 (x, y, z, 'Color', color_o, 'LineWidth', 1.5, 'HandleVisibility', 'off');

    rr_sat_old = rr_sat;

    % get frame for the video
    frame = getframe(gcf);
    writeVideo (video, frame);

    % set satellite orbital elements for the next iteration
    index = j + find (k == man_point);

    if ~isempty (index)
        a_s =       a(index);
        e_s =       e(index);
        i_s =       i(index);
        OM_s =      OM(index);
        om_s =      om(index);
        color_o =   color(index);
    end
end

close(video)



function theta = GetTrueAnomaly (thi, thf, dt, a, e, mu)
    % 
    % This function returns true anomaly values starting from 
    % 'thi' (initial true anomaly) to 'thf' (final true anomaly).
    % Those values are obtained by solving the following
    % implicit equation through Newton's method: 
    %       f(E) = t - sqrt( a^3/mu ) * (E - e*sin(E)) = 0
    %

    ti = TOF(a, e, 0, thi, mu); % initial time
    tf = TOF(a, e, 0, thf, mu); % final time
    T = 2*pi*sqrt( a^3/mu );    % orbit period

    while tf <= ti 
        tf = tf + T;
    end    

    time = ti+dt : dt : tf-dt;
    len = length (time);

    theta = zeros(1, len);
    E_old = 2*atan( sqrt( (1-e)/(1+e) ) * tan(thi/2));  % initial eccentric anomaly

    for k = 1:len 
        fun = @(E) time(k) - sqrt( a^3/mu ) * (E - e*sin(E));
        E_new = fzero(fun, E_old);  % solve the equation through Newton's method to obtain the eccentric anomaly

        theta(k) = 2*atan( sqrt( (1+e)/(1-e) ) * tan(E_new/2) );    % get true anomaly value
        E_old = E_new;
    end

    theta = [thi, theta, thf];
    theta = correctPeriod(theta);

end
