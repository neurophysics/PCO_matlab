function[vlen, vlen_der] = PCOa_obj_der(w, a, y, sign)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% start calculation of mean vector length and its partial derivatives     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filt = w'*y;
phase = angle(filt);

%% Result of the objective function #
a_norm = (a - mean(a))/std(a, 1);
sum1_no_square = mean(a_norm.*cos(phase));
sum2_no_square = mean(a_norm.*sin(phase));
% multiply with sign, i.e. -1 if function should be minimized or 1
% if it should be maximized
vlen = sign*sqrt(sum1_no_square^2 + sum2_no_square^2);

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
