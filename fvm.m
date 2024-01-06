T = 20;
delta_t = 0.01;
v = [1,2];
mu = 1/10;
u_bound = 0;
u_t0 = @(x,y) sin(pi*(x^2+y^2))*((x-1)^2+y^2-9);

[tr_nodes,voi_nodes,tr_elems,voi_elems] = mesh_generator(0.12);
new_tr_nodes = find_bound_nodes(tr_nodes,voi_nodes,voi_elems);
tr_col = size(new_tr_nodes,1);
U = zeros(tr_col,T);

for ui = 1:tr_col
    U(ui,1) = u_t0(tr_nodes(ui,1),tr_nodes(ui,2));
end

% area = calculate_area(50,voi_nodes,voi_elems)
% edge_length = calculate_edge_length(50, 71, voi_nodes, voi_elems)
% neighbour_matrix = find_neighbour(50, new_tr_nodes, voi_elems)
% node_length = calculate_node_length(50,51,new_tr_nodes)
% v_ij = calculate_vij (v,50,51,new_tr_nodes)

for timecount = 1:T
    for ui = 1:tr_col
        if new_tr_nodes(ui,3)==1
            U(ui,timecount+1) = u_bound;
        else
            neighbour_matrix = find_neighbour(ui, new_tr_nodes, voi_elems);
            sum_int = 0;
            sum_deter_Ui=0;
            sum_deter_Uj=0;
            for ni = 1: size(neighbour_matrix,1)
                if neighbour_matrix(ni,1) ~= 0 
                    % calculate inportant indexes
                    Ki = calculate_area(ui,voi_nodes,voi_elems);
                    Ki_inter_Kj = calculate_edge_length(ui,neighbour_matrix(ni,1),voi_nodes,voi_elems);
                    Xj_minus_Xi = calculate_node_length(ui,neighbour_matrix(ni,1),new_tr_nodes);
                    [v_ij_plus,v_ij_minus] = calculate_vij(v,ui,neighbour_matrix(ni,1),new_tr_nodes);

                    sum_deter_Ui = sum_deter_Ui - delta_t / Ki * Ki_inter_Kj * (v_ij_plus +mu/Xj_minus_Xi) *U(ui,timecount);
                    sum_deter_Uj = sum_deter_Uj - delta_t / Ki * Ki_inter_Kj * (v_ij_minus-mu/Xj_minus_Xi) *U(neighbour_matrix(ni,1),timecount);
                end
            end
            U(ui,timecount+1) = U(ui,timecount) + sum_deter_Ui + sum_deter_Uj + sum_int;
        end
    end
    time = delta_t * timecount;
    fprintf('%s %0.2f\n', "time =", time);
    % show voi plot
    figure;
    trisurf(tr_elems, tr_nodes(:,1), tr_nodes(:,2), U(:,timecount+1));
    zlim([-8, 8]);
    colorbar;
    clim([-8, 8]);
    drawnow;
    filename = sprintf('time_plot_%0.2f.png', time);
    %saveas(gcf, filename);
end

%==========================================================================
function [v_ij_plus,v_ij_minus] = calculate_vij (v,i,j,new_tr_nodes)
v_ij_plus = 0;
v_ij_minus = 0;
v_ij = v * [new_tr_nodes(j,1)-new_tr_nodes(i,1);
            new_tr_nodes(j,2)-new_tr_nodes(i,2)] /...
           calculate_node_length(i,j,new_tr_nodes);
if v_ij >= 0
    v_ij_plus = v_ij;
else
    v_ij_minus = v_ij;
end
end

%==========================================================================
function node_length = calculate_node_length(m,n,new_tr_nodes)
   node_length = sqrt((new_tr_nodes(m,1)-new_tr_nodes(n,1))^2+...
                      (new_tr_nodes(m,2)-new_tr_nodes(n,2))^2); 
end

%==========================================================================
function neighbour_matrix = find_neighbour(m, new_tr_nodes, voi_elems)
if new_tr_nodes(m,3) == 1
    neighbour_matrix =  ("This is a boundary point!");
else   
    c = size(new_tr_nodes,1);
    matrix = zeros(12);
    count = 0;

    for i = 1:c
        if i ~= m
            colm = size(voi_elems{m},2);
            coln = size(voi_elems{i},2);

            for j = 1:colm
                for k = 1:coln
                    if voi_elems{m}(j) == voi_elems{i}(k)
                       count = count+1;
                       matrix(count) = i;
                    end
                end
            end
        end
    end
    neighbour_matrix = unique(matrix);
end
end

%==========================================================================
function edge_length = calculate_edge_length(m, n, voi_nodes, voi_elems)
    colm = size(voi_elems{m},2);
    coln = size(voi_elems{n},2);
    elem_length_matrix = zeros(2,2);
    count = 0;
    for i = 1:colm
        for j = 1:coln
            if voi_elems{m}(i) == voi_elems{n}(j)
                count = count + 1;
                elem_length_matrix(count,1) = voi_nodes(voi_elems{m}(i),1);
                elem_length_matrix(count,2) = voi_nodes(voi_elems{m}(i),2);
            end
        end
    end
    edge_length = sqrt((elem_length_matrix(2,1) - elem_length_matrix(1,1))^ 2 + ...
                       (elem_length_matrix(2,2) - elem_length_matrix(1,2))^ 2);

    % m_matrix = zeros(colm,2);
    % for i = 1:colm
    %     m_matrix(i,1) = voi_nodes(voi_elems{m}(i),1);
    %     m_matrix(i,2) = voi_nodes(voi_elems{m}(i),2);
    % end
    % 
    % coln = size(voi_elems{n},2);
    % n_matrix = zeros(coln,2);
    % for i = 1:coln
    %     n_matrix(i,1) = voi_nodes(voi_elems{n}(i),1);
    %     n_matrix(i,2) = voi_nodes(voi_elems{n}(i),2);
    % end
    
    if edge_length == 0
        edge_length = ("There is no edge between these two points!");
    end
end

%==========================================================================
function area = calculate_area(m,voi_nodes,voi_elems)
    col = size(voi_elems{m},2);
    point_matrix = zeros(col,2);
    for i = 1:col
        point_matrix(i,1) = voi_nodes(voi_elems{m}(i),1);
        point_matrix(i,2) = voi_nodes(voi_elems{m}(i),2);
    end
    [~,area] = convhulln(point_matrix);
end

function [tr_nodes,voi_nodes,tr_elems,voi_elems] = mesh_generator(size)
    t = pi/12:pi/12:2*pi;
    pgon = polyshape({cos(t)*3+1, cos(t)}, ...
        {sin(t)*3, sin(t)});
    tr = triangulation(pgon);
    model = createpde;
    tnodes = tr.Points';
    telements = tr.ConnectivityList';
    
    geometryFromMesh(model,tnodes,telements);
    
    mesh = generateMesh(model,"Hmax", ...
    size, ...
           "GeometricOrder","linear");
    
    % show the triangle PDE mesh
    figure;
    pdemesh(mesh);

    tr_nodes = mesh.Nodes';
    tr_elems = mesh.Elements';
    [voi_nodes,voi_elems] = voronoin(tr_nodes);

    %show voi plot
    voronoi(tr_nodes(:,1),tr_nodes(:,2))
    axis equal
end

%==========================================================================
function new_tr_nodes = find_bound_nodes(tr_nodes,voi_nodes,voi_elems)
    row = size(tr_nodes,1);
    bd_matrix = zeros(row,1);
    for i =1:row
        col = size(voi_elems{i},2);
        if col == 4
            bd_matrix(i) = 1;
        else
            for j = 1:col
                if voi_nodes(voi_elems{i}(j),1)^2 + ...
                    voi_nodes(voi_elems{i}(j),2)^2 < 0.5^2
                    bd_matrix(i) = 1;
                elseif (voi_nodes(voi_elems{i}(j),1)-1)^2 + ...
                    voi_nodes(voi_elems{i}(j),2)^2 > 3^2
                    bd_matrix(i) = 1;
                end
            end
        end
    end
    new_tr_nodes = [tr_nodes, bd_matrix];
end