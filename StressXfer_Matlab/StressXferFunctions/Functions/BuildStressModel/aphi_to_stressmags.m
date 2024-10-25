% Calculate stress magnitudes from A_phi assuming frictional equilibrium
% Jens-Erik Lund Snee
% 2022

function [SHmax, Shmin, Sv, Pp] = aphi_to_stressmags(A_phi, z_km, mu)

% INPUTS:
z = z_km * 1000; % depth in m from input in km
g = 9.8; % m/s/s
rho = 2.6; % g/cm^3 rock density
rho = rho * (100^3)/1000; % kg/m^3 rock density
Sv = rho * g * z / 1000000; % MPa vertical stress at selected depth 
rho_water = 1; % g/cm^3 water density
rho_water = rho_water * (100^3)/1000; % kg/m^3 water density
Pp = rho_water * g * z / 1000000; % MPa hydrostatic pore pressure at selected depth
C0 = 0; % cohesion
mu_constraint = ((mu^2 + 1)^(1/2) + mu)^2;
%Sv = 26.01368*z/1000; % MPa/km value for 1.15 psi/ft


% convert A_phi to principal stresses assuming critical stress, hydrostatic, and constant overburden:
n = ones(length(A_phi),1); % preallocate n vector
phi = ones(length(A_phi),1); % preallocate phi vector

for i = 1:length(A_phi) % sweep over entire A_phi (column) dataset
    
    % find n: 
    if A_phi(i) <= 1
        n(i) = 0;
    elseif A_phi(i) >= 2
        n(i) = 2;
    else
        n(i) = 1;
    end
    
    phi(i) = (A_phi(i) - n(i) - 0.5) / ((-1)^n(i)) + 0.5; % find phi for each A_phi
    
    % find S1, S2, S3 for each A_phi value:
    if n(i) == 0 % normal faulting: Sv = S1, SHmax = S2, Shmin = S3:
        Shmin_vector(i) = (Sv-Pp) / mu_constraint + Pp;
        SHmax_vector(i) = (Sv-Shmin_vector(i)) * phi(i) + Shmin_vector(i);
    elseif n(i) == 1 % strike-slip faulting: SHmax = S1, Sv = S2, Shmin = S3:
        SHmax_vector(i) = ...
            ( (1/phi(i))*Sv - (1/phi(i))*Pp + Pp + Pp*((1/(mu_constraint*phi(i)))...
            - (1/mu_constraint)) ) / ( 1 + (1/(mu_constraint*phi(i)) - 1/mu_constraint) );
        Shmin_vector(i) = (SHmax_vector(i)-Pp)/mu_constraint + Pp;
    else % reverse faulting: SHmax = S1, Shmin = S2, Sv = S3:
        SHmax_vector(i) = mu_constraint * (Sv-Pp) + Pp;
        Shmin_vector(i) = phi(i) * (SHmax_vector(i)-Sv) + Sv;
    end
 
    % fill a stress tensor with the principal stresses for this iteration:
    principal_tensor{i} = [SHmax_vector(i) 0 0;
         0 Shmin_vector(i) 0;
         0 0 Sv];
    
end 

SHmax = SHmax_vector';
Shmin = Shmin_vector'; 

end



function [SHmax, Shmin] = aphi_to_hStresses(A_phi, mu, Sv, Pp)

% INPUTS:
%z = z_km * 1000; % depth in m from input in km
%g = 9.8; % m/s/s
%rho = 2.6; % g/cm^3 rock density
%rho = rho * (100^3)/1000; % kg/m^3 rock density
%Sv = rho * g * z / 1000000; % MPa vertical stress at selected depth 
%rho_water = 1; % g/cm^3 water density
%rho_water = rho_water * (100^3)/1000; % kg/m^3 water density
%Pp = rho_water * g * z / 1000000; % MPa hydrostatic pore pressure at selected depth
%C0 = 0; % cohesion
mu_constraint = ((mu^2 + 1)^(1/2) + mu)^2;
%Sv = 26.01368*z/1000; % MPa/km value for 1.15 psi/ft


% convert A_phi to principal stresses assuming critical stress, hydrostatic, and constant overburden:
n = ones(length(A_phi),1); % preallocate n vector
phi = ones(length(A_phi),1); % preallocate phi vector

for i = 1:length(A_phi) % sweep over entire A_phi (column) dataset
    
    % find n: 
    if A_phi(i) <= 1
        n(i) = 0;
    elseif A_phi(i) >= 2
        n(i) = 2;
    else
        n(i) = 1;
    end
    
    phi(i) = (A_phi(i) - n(i) - 0.5) / ((-1)^n(i)) + 0.5; % find phi for each A_phi
    
    % find S1, S2, S3 for each A_phi value:
    if n(i) == 0 % normal faulting: Sv = S1, SHmax = S2, Shmin = S3:
        Shmin_vector(i) = (Sv-Pp) / mu_constraint + Pp;
        SHmax_vector(i) = (Sv-Shmin_vector(i)) * phi(i) + Shmin_vector(i);
    elseif n(i) == 1 % strike-slip faulting: SHmax = S1, Sv = S2, Shmin = S3:
        SHmax_vector(i) = ...
            ( (1/phi(i))*Sv - (1/phi(i))*Pp + Pp + Pp*((1/(mu_constraint*phi(i)))...
            - (1/mu_constraint)) ) / ( 1 + (1/(mu_constraint*phi(i)) - 1/mu_constraint) );
        Shmin_vector(i) = (SHmax_vector(i)-Pp)/mu_constraint + Pp;
    else % reverse faulting: SHmax = S1, Shmin = S2, Sv = S3:
        SHmax_vector(i) = mu_constraint * (Sv-Pp) + Pp;
        Shmin_vector(i) = phi(i) * (SHmax_vector(i)-Sv) + Sv;
    end
 
    % fill a stress tensor with the principal stresses for this iteration:
    principal_tensor{i} = [SHmax_vector(i) 0 0;
         0 Shmin_vector(i) 0;
         0 0 Sv];
    
end 

SHmax = SHmax_vector;
Shmin = Shmin_vector; 

end