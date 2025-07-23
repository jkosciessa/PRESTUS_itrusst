% Parameters
f = 1e6;                % Frequency [Hz]
c = 1480;               % Speed of sound in water [m/s]
lambda = c / f;         % Wavelength [m]
k = 2 * pi / lambda;    % Wavenumber
a = 10e-3;              % Piston radius [m]
N_r = 100;              % Number of radial divisions
N_theta = 200;          % Number of angular divisions

% Observation grid (axial cut)
z = linspace(0.01, 0.1, 300);    % Axial positions [m]
x = linspace(-0.02, 0.02, 200);  % Lateral positions [m]
[X, Z] = meshgrid(x, z);
Y = zeros(size(X));              % y=0 plane

% Discretize transducer surface (polar coordinates)
r = linspace(0, a, N_r);
theta = linspace(0, 2*pi, N_theta);
[RR, TT] = meshgrid(r, theta);
dA = (a/N_r) * (2*pi/N_theta) * RR; % Area element at each (r,theta)

% Transducer surface points
x0 = RR .* cos(TT);
y0 = RR .* sin(TT);

% Initialize pressure field
P = zeros(size(X));

% Loop over transducer surface elements
for idx_r = 1:N_r
    for idx_theta = 1:N_theta
        % Source point on transducer
        xs = x0(idx_theta, idx_r);
        ys = y0(idx_theta, idx_r);
        zs = 0;
        
        % Distance from source to field point
        R = sqrt((X - xs).^2 + (Y - ys).^2 + (Z - zs).^2);
        
        % Rayleigh-Sommerfeld integrand (monopole, rigid baffle)
        integrand = (Z ./ R) .* exp(1i * k * R) ./ R;
        
        % Sum contribution (midpoint rule)
        P = P + integrand * dA(idx_theta, idx_r);
    end
end

% Normalize (optional)
P = P / max(abs(P(:)));

% Plot intensity profile
imagesc(x*1e3, z*1e3, 20*log10(abs(P)));
xlabel('Lateral position x (mm)');
ylabel('Axial position z (mm)');
title('Acoustic Pressure Field (dB, Rayleigh-Sommerfeld)');
colorbar; colormap('jet');
axis xy;
