function [vlen, wy] = PCO_a (a, y, nu, bestof)

%define nested functions
    function [obj, obj_der] = objfun_y(x)
        [obj, obj_der] = PCOa_obj_der(x, a, y, -1);
    end

    function [obj, obj_der] = objfun_by(x)
        [obj, obj_der] = PCOa_obj_der(x, a, by, -1);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Phase Amplitude Coupling Optimization, variant with provided amplitude
%
%
%     It maximizes the length of the "mean vector" and
%     returns the filter coefficients w
%
%     Input:
%     ------
%     a - (1d numpy array, floats > 0) amplitudes
%     Y - (2d numpy array, complex) analytic representation of signal,
%         channels x datapoints
%     num - (int > 0) - determines the number of filters that will be
%                       derived. This depends also on the rank of Y,
%                       the final number of filters will be min
%                       ([num, rank(Y)]), defaults to 1
%
%     bestof (int > 0) - the number of restarts for the optimization of the
%                        individual filters. The best filter over all
%                        these restarts with random initializations is
%                        chosen, defaults to 15.
%
%     Output:
%     -------
%     vlen - numpy array - the length of the mean vector for each filter
%     Wy - numpy array - the filters for Y, each filter is in a column of
%     Wy (if num==1: Wx is 1d)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fill in unset optional values.
switch nargin
    case 2
        nu = 1;
        bestof = 15;
    case 3
        bestof = 15;
end
% define the options of the minimizer functions
minoptions = optimoptions('fminunc', ...
    'Algorithm', 'quasi-newton', ...
    'HessUpdate', 'bfgs', ...
    'GradObj','on', ...
    'Display', 'off');

%%Whiten the (real part of the) data
[why, sy, ~] = svd(real(y), 'econ');
%get rank
py = rank(y);

why = bsxfun(@rdivide, why(1:end,1:py)', diag(sy(1:py,1:py)))';
%whiten for the real part of the data
y = why'*y;

%get the final number of filters
num = min([nu, py]);

%% start optimization
for i = 1:num    
    if i == 1
        % get first filter
        % get best parameters and function values of each run
        all_best = 0;
        for k = 1:bestof
            [x_best, best_f] = fminunc(@objfun_y, rand(py,1)*2 - 1, ...
                minoptions);
            if best_f < all_best
                all_best = best_f;
                all_best_x = x_best;
            end
        end
        %save results
        vlen = all_best;
        filt = all_best_x;
    else
        %get consecutive pairs of filters
        %project data into null space of previous filters
        %this is done by getting the right eigenvectors of the filter
        %maxtrix corresponding to vanishing eigenvalues
        [~, ~, Vy] = svd(filt');
        by = Vy(1:end,i:end)'*y;
        % get best parameters and function values of each run
        all_best = 0;
        for k = 1:bestof
            [x_best, best_f] = fminunc(@objfun_by, rand(py - i + 1,1)*2 - 1, ...
                minoptions);
            if best_f < all_best
                all_best = best_f;
                all_best_x = x_best;
            end
        end
        % save results
        vlen = cat(2, vlen, all_best);
        filt = cat(2, filt, Vy(1:end,i:end)*all_best_x);
    end
    %project filters back into original (un-whitened) channel space
    wy = why*filt;
    %normalize filters to have unit length
    wy = bsxfun(@rdivide, wy, sqrt(sum(wy.^2, 1)));
    vlen= -1*vlen;
end
end
