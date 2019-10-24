function C = setupEliminationTemplates(data)
    [coeffs,coeffs_ind] = computeCoefficients(data);

    C_ind = [1,6,7,8,9,12,13,14,15,18,19,22,23,24,25,26,27,28,29,30,32,33,40,41,43,46,47,48,50,51,52,53,56,57,58,59];
    C = zeros(6,10);
    C(C_ind) = coeffs(coeffs_ind);
end
