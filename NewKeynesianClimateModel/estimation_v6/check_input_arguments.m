function check_input_arguments(options_, M_, oo_)
%function check_input_arguments(options_, M_, oo_)
%Conducts checks for inconsistent/missing inputs to deterministic
%simulations
% Inputs:
%   options_            [structure] describing the options
%   M_                  [structure] describing the model
%   oo_                 [structure] storing the results

% Copyright © 2015-2023 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <https://www.gnu.org/licenses/>.

if options_.stack_solve_algo < 0 || options_.stack_solve_algo > 7
    error('perfect_foresight_solver:ArgCheck','PERFECT_FORESIGHT_SOLVER: stack_solve_algo must be between 0 and 7')
end

if ~options_.block && ~options_.bytecode && options_.stack_solve_algo ~= 0 ...
        && options_.stack_solve_algo ~= 1 && options_.stack_solve_algo ~= 6 ...
        && options_.stack_solve_algo ~= 7
    error('perfect_foresight_solver:ArgCheck','PERFECT_FORESIGHT_SOLVER: you must use stack_solve_algo={0,1,6,7} when not using block nor bytecode option')
end

if options_.block && ~options_.bytecode && options_.stack_solve_algo == 5
    error('perfect_foresight_solver:ArgCheck','PERFECT_FORESIGHT_SOLVER: you can''t use stack_solve_algo = 5 without bytecode option')
end

if isempty(oo_.endo_simul) || any(size(oo_.endo_simul) ~= [ M_.endo_nbr, M_.maximum_lag+options_.periods+M_.maximum_lead ])

    if options_.initval_file
        fprintf('PERFECT_FORESIGHT_SOLVER: ''oo_.endo_simul'' has wrong size. Check whether your initval-file provides %d periods.',M_.maximum_endo_lag+options_.periods+M_.maximum_endo_lead)
        error('perfect_foresight_solver:ArgCheck','PERFECT_FORESIGHT_SOLVER: ''oo_.endo_simul'' has wrong size. Did you run ''perfect_foresight_setup'' ?')
    else
        error('perfect_foresight_solver:ArgCheck','PERFECT_FORESIGHT_SOLVER: ''oo_.endo_simul'' has wrong size. Did you run ''perfect_foresight_setup'' ?')
    end
end

if (M_.exo_nbr > 0) && ...
        (isempty(oo_.exo_simul) || any(size(oo_.exo_simul) ~= [ M_.maximum_lag+options_.periods+M_.maximum_lead, M_.exo_nbr ]))
    if options_.initval_file
        fprintf('PERFECT_FORESIGHT_SOLVER: ''oo_.exo_simul'' has wrong size. Check whether your initval-file provides %d periods.',M_.maximum_endo_lag+options_.periods+M_.maximum_endo_lead)
        error('perfect_foresight_solver:ArgCheck','PERFECT_FORESIGHT_SOLVER: ''oo_.exo_simul'' has wrong size.')
    else
        error('perfect_foresight_solver:ArgCheck','PERFECT_FORESIGHT_SOLVER: ''oo_.exo_simul'' has wrong size. Did you run ''perfect_foresight_setup'' ?')
    end
end
