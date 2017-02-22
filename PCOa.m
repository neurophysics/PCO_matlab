function [vlen, wy] = PCOa(a, y, nu, bestof)
    % Phase Amplitude Coupling Optimization, variant with provided
    % amplitude
    %
    % It maximizes the length of the "mean vector" and returns the filter
    % coefficients wy
    %
    % Input:
    % ------
    % a - (column vector) amplitudes
    % y - (2d array, complex) analytic representation of signal,
    %     channels x datapoints
    % num - (int > 0) - determines the number of filters that will be
    %                   derived. This depends also on the rank of y,
    %                   the final number of filters will be min
    %                   ([num, rank(y)]), defaults to 1
    % bestof (int > 0) - the number of restarts for the optimization of the
    %                    individual filters. The best filter over all
    %                    these restarts with random initializations is
    %                    chosen, defaults to 15.
    %
    % Output:
    % -------
    % vlen - row vector - the length of the mean vector for each filter
    % wy - 2d array - the filters for Y, each filter is in a column of Wy
    %define nested functions for the optimization
    %for the first filter
    function [obj, obj_der] = objfun_y(x)
        [obj, obj_der] = PCOa_obj_der(x, a, y, -1);
    end
    %for subsequent filters
    function [obj, obj_der] = objfun_by(x)
        [obj, obj_der] = PCOa_obj_der(x, a, by, -1);
    end

    % Fill in standar values for unset optional values.
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

    %Whiten the (real part of the) data
    [why, sy, ~] = svd(real(y), 'econ');
    %get rank
    py = rank(y);
    % calculate the whitening filter
    why = bsxfun(@rdivide, why(1:end,1:py)', diag(sy(1:py,1:py)))';
    % apply the whitening filter
    y = why'*y;

    %get the final number of filters
    num = min([nu, py]);

    % start optimization
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
                [x_best, best_f] = fminunc(@objfun_by, ...
                    rand(py - i + 1,1)*2 - 1, minoptions);
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

% The following function defines the objective function and its
% derivative
function[vlen, vlen_der] = PCOa_obj_der(w, a, y, sign)
    % Calculation of mean vector length and its partial derivatives
    % at the filter y
    % if sign == -1, the results are multiplied with -1, this enables
    % maximizing the original variable (by minimizing -1*the objective)
    filt = w'*y;
    phase = angle(filt);
    % calculate the result of the objective function
    a_norm = (a - mean(a))/std(a, 1);
    sum1_no_square = mean(a_norm.*cos(phase));
    sum2_no_square = mean(a_norm.*sin(phase));
    % multiply with sign, i.e. -1 if function should be minimized
    % or 1 if it should be maximized
    vlen = sign*sqrt(sum1_no_square^2 + sum2_no_square^2);

    % the gradient should only be calculated if nargout > 1 
    if nargout > 1 % calculate gradient
        % Partial derivative of phase
        phase_dwi = bsxfun(@rdivide, ...
            bsxfun(@times, -real(y), imag(filt)) + ...
            bsxfun(@times, imag(y), real(filt)), ...
            real(filt).^2 + imag(filt).^2);
        % Derivative of first summand
        sum1_d = 2*sum1_no_square*mean(bsxfun(@times, phase_dwi, ...
            -a_norm.*sin(phase)), 2);
        % Derivative of second summand
        sum2_d = 2*sum2_no_square*mean(bsxfun(@times, phase_dwi, ...
            a_norm.*cos(phase)), 2);
        % Derivative of sum
        sum_d = sum1_d + sum2_d;
        % Derivative including square root
        vlen_der = sum_d/(2*vlen);
    end
end

