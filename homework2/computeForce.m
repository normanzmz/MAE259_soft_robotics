function [Fe, Je] = computeForce(x, m1, m2, refTwist)

global EI EA GJ voronoiRefLen kappaBar refLen 

ndof = numel(x);
ne = (ndof+1)/4 - 1; % number of edges

Fe = zeros(ndof, 1);
Je = zeros(ndof, ndof);

for c = 2:ne
    ind = 4*c-7 : 4*c+3;
    node0 = [x(4*c-7), x(4*c-6), x(4*c-5)]; 
    node1 = [x(4*c-3), x(4*c-2), x(4*c-1)]; 
    node2 = [x(4*c+1), x(4*c+2), x(4*c+3)];
    theta_e = x(4*c-4);
    theta_f = x(4*c);
    l_k = voronoiRefLen(c);    
    m1e = m1(c-1,:);
    m2e = m2(c-1,:);
    m1f = m1(c,:);
    m2f = m2(c,:);

    [dF, dJ] = gradEb_hessEb(node0, node1, node2, ...
    m1e, m2e, m1f, m2f, kappaBar(c,:), l_k, EI);    
    Fe(ind) = Fe(ind) - dF;
    Je(ind, ind) = Je(ind, ind) - dJ;

    [dF, dJ] = gradEt_hessEt(node0, node1, node2, ...
        theta_e, theta_f, refTwist(c), l_k, GJ);
    Fe(ind) = Fe(ind) - dF;
    Je(ind, ind) = Je(ind, ind) - dJ;
end

for c = 1:ne
    ind = 4*c-3 : 4*c+3;
    node0 = [x(4*c-3), x(4*c-2), x(4*c-1)]; 
    node1 = [x(4*c+1), x(4*c+2), x(4*c+3)]; 
    l_k = refLen(c);    
    [dF, dJ] = gradEs_hessEs(node0, node1, l_k, EA);
    Fe(ind) = Fe(ind) - dF;
    Je(ind, ind) = Je(ind, ind) - dJ;
end

% end