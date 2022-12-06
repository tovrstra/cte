% All values are in SI base units.
global STD_GRAVITY
global COULOMB_CONSTANT
global MASS
global LENGTH
global MOMENT_OF_INERTIA
global TIMESTEP
global NSTEP

STD_GRAVITY = 9.80665;
COULOMB_CONSTANT = 8.9875517923e9;
MASS = 0.062;
LENGTH = 0.80;
MOMENT_OF_INERTIA = MASS * LENGTH^2;
TIMESTEP = 1.79 / 50;
NSTEP = 75;

% Initialize some variables
theta = pi/3;
omega = 0.0;

% Initialize zero arrays in which the trajectory is stored.
thetas = 1: NSTEP;
omegas = 1: NSTEP;

% Store the initial state
thetas(1) = theta;
omegas(1) = omega;

% Leapfrog integration
ang_accel = torque_theta(theta) / MOMENT_OF_INERTIA;
for istep = 1: NSTEP
    % Leapfrog algorithm for a single step.
    omega_half = omega + TIMESTEP * ang_accel / 2;
    theta_next = theta + TIMESTEP * omega_half;
    ang_accel_next = torque_theta(theta_next) / MOMENT_OF_INERTIA;
    omega_next = omega_half + TIMESTEP * ang_accel_next / 2;

    % Store results.
    thetas(istep + 1) = theta_next;
    omegas(istep + 1) = omega_next;

    % Assign "next" results to current.
    theta = theta_next;
    omega = omega_next;
    ang_accel = ang_accel_next;
end

% Plot of position, velocity and energies as a function of time.

ts = (1 : NSTEP + 1) * TIMESTEP;

hold on;
subplot(1, 3, 1)
plot(ts, thetas)
xlabel("time t [s]")
ylabel("angle theta [rad]")

subplot(1, 3, 2)
plot(ts, omegas)
xlabel("time t [s]")
ylabel("angular velocity omega [rad/s]")

subplot(1, 3, 3)
eks = kinetic_energy(omegas);
eps = potential_energy(thetas);
plot(ts, eks, "b", ts, eps, "r", ts, eks + eps, "k")
% plot(ts, eps)
% plot(ts, eks + eps, 'color', "k")
legend("kin. E", "pot. E", "total E")
xlabel("time t [s]")
ylabel("energy [J]")
saveas(gca(), "matlab_trajectory.png", "png")
clf()
hold off;


% Plot velocity versus position

hold on;
plot(thetas, omegas)
xlabel("angle theta [rad]")
ylabel("angular velocity omega [rad/s]")
saveas(gca(), "matlab_omega_theta.png", "png")
clf()
hold off;


% Plot contour lines for the total energy

theta_low = min(thetas);
theta_high = max(thetas);
theta_margin = (theta_high - theta_low) * 0.15;
theta_low = theta_low - theta_margin;
theta_high = theta_high + theta_margin;
omega_low = min(omegas);
omega_high = max(omegas);
omega_margin = (omega_high - omega_low) * 0.15;
omega_low = omega_low - omega_margin;
omega_high = omega_high + omega_margin;
hold on;
[theta_grid, omega_grid] = meshgrid(linspace(theta_low, theta_high, 51), linspace(omega_low, omega_high, 51));
energy_grid = kinetic_energy(omega_grid) + potential_energy(theta_grid);
contour(theta_grid, omega_grid, energy_grid, 40)
xlabel("angle theta [rad]")
ylabel("angular velocity omega [rad/s]")
saveas(gca(), "matlab_ener_contour.png", "png")
clf()
hold off;

% Functions must be defined in the end, in Matlab,
% unlike in any other programming language...

function [out] = torque_theta(theta)
    global MASS
    global STD_GRAVITY
    global LENGTH
    out = -MASS * STD_GRAVITY * LENGTH * sin(theta);
end

function [out] = potential_energy(theta)
    global MASS
    global STD_GRAVITY
    global LENGTH
    out = MASS * STD_GRAVITY * LENGTH * (1 - cos(theta));
end

function [out] = kinetic_energy(omega)
    global MOMENT_OF_INERTIA
    out = (0.5 * MOMENT_OF_INERTIA) * omega.^2;
end
