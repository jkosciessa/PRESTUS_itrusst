% Debug script for k-Wave grid size in mm
% Assume kgrid is your kWaveGrid object already created

% Display all relevant grid properties
fprintf('Grid points: Nx=%d, Ny=%d, Nz=%d\n', kgrid.Nx, kgrid.Ny, kgrid.Nz);
fprintf('Grid spacing (m): dx=%.3e, dy=%.3e, dz=%.3e\n', kgrid.dx, kgrid.dy, kgrid.dz);
fprintf('Grid spacing (mm): dx=%.3f, dy=%.3f, dz=%.3f\n', kgrid.dx*1e3, kgrid.dy*1e3, kgrid.dz*1e3);

% Physical sizes
x_size_mm = kgrid.x_size * 1e3;
y_size_mm = kgrid.y_size * 1e3;
z_size_mm = kgrid.z_size * 1e3;
fprintf('Physical sizes (mm): x=%.3f, y=%.3f, z=%.3f\n', x_size_mm, y_size_mm, z_size_mm);

% Total grid points
fprintf('Total grid points: %d\n', kgrid.total_grid_points);

% Verify calculations
fprintf('Verification - x_size calc: %d * %.3e m = %.3f mm (matches: %s)\n', ...
    kgrid.Nx, kgrid.dx, kgrid.Nx * kgrid.dx * 1e3, ...
    string(abs(kgrid.Nx * kgrid.dx * 1e3 - x_size_mm) < 1e-10));

% Inspect full object properties (for deeper debugging)
disp(kgrid);  % Shows all properties
