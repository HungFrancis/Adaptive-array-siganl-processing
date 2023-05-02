function calc_nb_beampattern(b)
%CALC_NB_BEAMPATTERN(b) Calculate the beamformer coefficients
%
% Calculate the beampattern coefficient for each angle that is stored in
% b.angles.
%
% The narrowband beamformer weights are available through b.nb_weights
%
% If you need array response vectors you can calculate them using the
% function: b.array_response_vector(theta, f)
%
% Store the resulting beampattern in b.nb_beampattern.

b.nb_beampattern = abs((b.nb_weights'*b.array_response_vector(b.angles, b.nb_frequency)).^2)/(b.array.number_of_sensors^2);
end

