% author: A. Panichella
% title: An Adaptive Evolutionary Algorithm based on Non-Euclidean Geometry for Many-objective Optimization
% venue: GECCO 2019

function d = Point2LineDistance(P, A, B)
    d = zeros(size(P,1),1);
    for i = 1:size(P,1)
        pa = P(i,:) - A;
        ba = B - A;
        t = dot(pa, ba)/dot(ba, ba);
        d(i,1) = norm(pa - t * ba,2);
    end
end