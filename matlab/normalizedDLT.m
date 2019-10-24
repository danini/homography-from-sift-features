
%% The normalized DLT algorithm
function H = normalizedDLT(pts1, pts2)
    [normalized_src_points, normalized_dst_points, T1, T2] = normalizePoints(pts1, pts2);
    H = DLT(normalized_src_points, normalized_dst_points);
    H = inv(T2) * H * T1;
end