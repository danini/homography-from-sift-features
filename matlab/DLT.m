%% The DLT algorithm
function H = DLT(pts1,pts2)
    NUMP=size(pts1,1);

    A=[];
    for i=1:NUMP
        x1=pts1(i,1);
        y1=pts1(i,2);

        x2=pts2(i,1);
        y2=pts2(i,2);

        A(2*i-1,:)=[-x1,-y1,-1.0,0,0,0,x2*x1,x2*y1,x2];
        A(2*i,:)=[0,0,0,-x1,-y1,-1.0,y2*x1,y2*y1,y2];
    end

    [tmp1,tmp2,V]=svd(A);
    h = V(:,9);
    H = [h(1),h(2),h(3);h(4),h(5),h(6);h(7),h(8),h(9)];
end