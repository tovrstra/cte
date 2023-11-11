% All values are in SI base units.
global MASS
global FORCE_CONST
global TIMESTEP
global NSTEP

MASS = 0.062;
FORCE_CONST = 2.5;
TIMESTEP = 0.01;
NSTEP = 75;

% Initialize some variables
x = 0.3;
v = 0.0;

% Initialize zero arrays in which the trajectory is stored.
xs = 1: NSTEP;
vs = 1: NSTEP;

% Store the initial state
xs(1) = x;
vs(1) = v;

% Leapfrog integration
a = force(x) / MASS;
for istep = 1: NSTEP
    % Leapfrog algorithm for a single step.
    v_half = v + TIMESTEP * a / 2;
    x_next = x + TIMESTEP * v_half;
    a_next = force(x_next) / MASS;
    v_next = v_half + TIMESTEP * a_next / 2;

    % Store results.
    xs(istep + 1) = x_next;
    vs(istep + 1) = v_next;

    % Assign "next" results to current.
    x = x_next;
    v = v_next;
    a = a_next;
end

% Plot of position, velocity and energies as a function of time.

ts = (1 : NSTEP + 1) * TIMESTEP;

hold on;
subplot(1, 3, 1)
plot(ts, xs)
xlabel("t [s]")
ylabel("x [m]")

subplot(1, 3, 2)
plot(ts, vs)
xlabel("t [s]")
ylabel("v [m/s]")

subplot(1, 3, 3)
eks = kinetic_energy(vs);
eps = potential_energy(xs);
plot(ts, eks, "b", ts, eps, "r", ts, eks + eps, "k")
legend("kin. E", "pot. E", "total E")
xlabel("t [s]")
ylabel("Energy [J]")
saveas(gca(), "matlab_trajectory.png", "png")
clf()
hold off;


% Plot velocity versus position

hold on;
plot(xs, vs)
xlabel("x [m]")
ylabel("v [m/s]")
saveas(gca(), "matlab_v_x.png", "png")
clf()
hold off;


% Plot contour lines for the total energy

x_low = min(thetas);
x_high = max(thetas);
x_margin = (x_high - x_low) * 0.15;
x_low = x_low - x_margin;
x_high = x_high + x_margin;
v_low = min(omegas);
v_high = max(omegas);
v_margin = (v_high - v_low) * 0.15;
v_low = v_low - v_margin;
v_high = v_high + v_margin;
hold on;
[x_grid, v_grid] = meshgrid(linspace(x_low, x_high, 51), linspace(v_low, v_high, 51));
energy_grid = kinetic_energy(v_grid) + potential_energy(x_grid);
contour(x_grid, v_grid, energy_grid, 40)
xlabel("x [m]")
ylabel("v [m/s]")
saveas(gca(), "matlab_ener_contour.png", "png")
clf()
hold off;

% Functions must be defined in the end, in Matlab,
% unlike in any other programming language...

function [out] = force(x)
    global FORCE_CONST
    out = -FORCE_CONST * x;
end

function [out] = potential_energy(x)
    global FORCE_CONST
    out = (0.5 * FORCE_CONST) * x.^2;
end

function [out] = kinetic_energy(v)
    global MASS
    out = (0.5 * MASS) * v.^2;
end
