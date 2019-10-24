function [normalized_src_points, normalized_dst_points, T1, T2] = normalizePoints(pts1, pts2)
    N			= size(pts1,1);
    
	massPoint1	= mean(pts1);
	massPoint2	= mean(pts2);
    
    normalized_src_points  = pts1 - repmat(massPoint1, N, 1);
    normalized_dst_points  = pts2 - repmat(massPoint2, N, 1);
    
    avgDist1    = 0;
    avgDist2    = 0;
    
	for i = 1 : N
		avgDist1		= avgDist1 + norm(normalized_src_points(i,:));
		avgDist2		= avgDist2 + norm(normalized_dst_points(i,:));
    end
    
	avgDist1	= avgDist1 / N;
	avgDist2	= avgDist2 / N;
	avgRatio1	= sqrt(2) / avgDist1;
	avgRatio2	= sqrt(2) / avgDist2;
	
    normalized_src_points  = normalized_src_points * avgRatio1;
    normalized_dst_points  = normalized_dst_points * avgRatio2;
    normalized_src_points(:,3)	= 1;
    normalized_dst_points(:,3)	= 1;
    		
	T1			= [avgRatio1, 0, 0;
					0, avgRatio1, 0;
					0, 0, 1] * [1, 0, -massPoint1(1);
					0, 1, -massPoint1(2);
					0, 0, 1];
	
	T2			= [avgRatio2, 0, 0;
					0, avgRatio2, 0;
					0, 0, 1] * [1, 0, -massPoint2(1);
					0, 1, -massPoint2(2);
					0, 0, 1];
end

