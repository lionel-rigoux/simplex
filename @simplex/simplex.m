classdef simplex < handle
    
    properties (SetAccess = private)
        resolution = 20;
        %end
        
        % properties (SetAccess = private, GetAccess = private, Hidden = true)
        unit;
        grid;
        centers;
        facets;
        vertices;
    end
    
    %% ====================================================================
    % Static methods
    %  ====================================================================
        methods (Access = private, Hidden = true)
            function q =  normalize (s, q)
                q = bsxfun (@rdivide, q, sum (q, 2));
            end
        end
        
        methods (Access = public, Hidden = false)
        function xy = barycentricToCartesian (s, q)
            n = size (q, 1);
            xy = s.unit.barycentricToCartesian (ones(n,1), q);
        end
        
        function q = cartesianToBarycentric (s, xy)
            n = size (xy, 1);
            q = s.unit.cartesianToBarycentric (ones(n, 1), xy);
            q = max (q, 0);
            q = min (q, 1);
            q = s.normalize(q);
        end
    end
    
    methods (Access = private, Hidden = true)
        
        
        
        function [t, vertices] = getGrid (s, n)
            %steps = linspace (0, 1, n);
            
            vertices = [];
            for q1 = 0 : n
                for q2 = 0 : (n-q1)
                    q3 = n-q1-q2;
                    vertices(end+1,:) = [q1 q2 q3]/(q1+q2+q3);
                end
            end
            points = s.barycentricToCartesian (vertices);
            t = delaunayTriangulation(points);
            
        end
    end
    
    
    methods
        
        
        function s = simplex ()
            
            % reference frame
            s.unit = delaunayTriangulation([0 0; 1 0; 0.5 sin(pi/3)]);
            
            % grid
            
            [s.grid,s.vertices] = s.getGrid(s.resolution);
            
            % centers
            s.centers = s.cartesianToBarycentric(incenter (s.grid));
            
        end
        
        function p = plot_axis (s)
            hold on
            t(1) = patch('XData',s.unit.Points(:,1),'YData',s.unit.Points(:,2));
            t(1).FaceColor = [235 235 235]/255;
            t(1).EdgeColor = 'none';
            axis square
            axis off
            %t(1).ZData = -.01*ones(size(t(1).XData));
            
            
            t(2:3)=plot_grid (s);
            t(3).LineWidth = 2;
            t(3).Color = [1 1 1 .5];
            t(2).LineWidth = 1;
            t(2).Color = [1 1 1 .6];
            
            %t(2).ZData = -.01*ones(size(t(2).XData));
            %t(3).ZData = -.01*ones(size(t(3).XData));
            
        end
        
        function t=plot_grid (s)
            hold on
            min = getGrid (s, 20);
            t(1) = triplot (min, 'Color', .8*[1 1 1], 'LineWidth', .3);
            t(1).Color(4) = .3;
            maj = getGrid (s, 4);
            t(2)=triplot (maj, 'Color', .95*[1 1 1], 'LineWidth', .3);
            t(2).Color(4) = .4;
            
            %t(1).ZData = +eps*ones(size(t(1).XData));
            %t(2).ZData = +eps*ones(size(t(2).XData));
            
        end
        
        function h = plot_function (s, fHandle)
            
            
            val = fHandle(s.centers);
            
            for i = 1 : s.grid.size(1)
                pointsId = s.grid.ConnectivityList(i, :);
                coor.x(:, i) = s.grid.Points(pointsId, 1);
                coor.y(:, i) = s.grid.Points(pointsId, 2);
            end
            
            h = patch (coor.x, coor.y, val);
            h.EdgeColor = 'none';
            
            colormap (viridis ());
        end
        
        function h = scatter (s, qs, varargin)
            
            if nargin < 3
                varargin = {'Marker','o','MarkerFaceColor',[.8 0 0],'MarkerEdgeColor','none','SizeData',35};
            end
            
            points = s.barycentricToCartesian (qs);
            
            h = scatter(points(:,1), points(:,2), varargin{:});
            
            
        end
        
        function h = fcontour (s, fHandle, color)
            
            function v = interpolant(x,y)
                dims = size(x);
                q = s.cartesianToBarycentric([x(:),y(:)]);
                v = fHandle(q);
                % prevent extrapolation
                isIn = inpolygon (x(:), y(:), s.unit.Points(:,1), s.unit.Points(:,2));
                v(~isIn) = NaN;
                % avoid border border evaluation
                minQ = 1e-8;
                v(any(q < minQ, 2)) = NaN;
                v = reshape(v,dims);
            end
            
            % plot contour
            boundaries = [0 1 0 max(s.grid.Points(:,2))];
            h= fcontour(@interpolant, boundaries,'MeshDensity',300);
            
            % make pretty
            h.LineWidth=1.5;
            if exist('color','var')
                h.LineColor = color;
            end
            
            %clabel(C,h,'LabelSpacing',Inf); % ,'Color',h.LineColor,
            
            colormap (viridis ());
            
        end
        
        function h = patch (s, fHandle)
            
            points = s.barycentricToCartesian (s.vertices);
            x = points(:,1);
            y = points(:,2);
            h=patch('Faces',delaunay(x,y),'Vertices',[x y]);
            h.FaceVertexCData = fHandle(s.vertices);
            h.FaceColor = 'interp';
            h.EdgeColor = 'none';
            
            colormap (viridis ());
            
        end
        function h = surf (s, fHandle)
            
            points = s.barycentricToCartesian (s.vertices);
            x = points(:,1);
            y = points(:,2);
            
            v = fHandle(s.vertices);
            tri = delaunay(s.vertices(:,1:2));
            
            h = trimesh(tri,x,y,v);
            
            % evaluate function
            
            h.FaceColor = 'interp';
            
            h.EdgeColor = 'w';
            h.EdgeAlpha = 0.2;
            
            colormap (viridis ());
            
            
        end
        
        function [N, edges] = histcounts (s, qs, nBins, varargin)
            
            if nargin < 3
                N = histcounts2 (qs(:,1), qs(:,2));
                nBins = round(mean(size(N)));
            end
            
            edges = linspace (0, 1 ,nBins + 1);
            N = histcounts2 (qs(:,1), qs(:,2), edges, edges);
            
            isOut = fliplr(tril(ones(nBins)) - eye(nBins));
            N(isOut==1) = NaN;
            
            
        end
        
        function h = histogram (s, qs)
            
            points = s.barycentricToCartesian(qs);
            ids = s.grid.pointLocation (points);
            val = histcounts (ids, (1:s.grid.size(1)+1) - 0.5);
            %val = val / sum(val);
            
            points = s.barycentricToCartesian (s.vertices);
            x = points(:,1);
            y = points(:,2);
            h=patch('Faces',delaunay(x,y),'Vertices',[x y]);
            
            h.FaceColor = 'flat';
            h.EdgeColor = 'none';
            h.FaceVertexCData = val';
            
            
            
            colormap (viridis ());
            
            
        end
        
        function h = contourHist (s, qs, varargin)
            

            [N,c] = s.histcounts(qs(:,1:2),varargin{:});
            [X,Y] = meshgrid((c(1:end-1) + c(2:end))/2);
            
            F = griddedInterpolant(X',Y',N);
            
            function v = inter(q)
                v = F(q(:,1), q(:,2));
            end
            
            h=s.fcontour(@inter);
            
        end
        
        function h = quiver (s, q, vector)
            origin = s.barycentricToCartesian(q);
            endpoint = s.barycentricToCartesian(s.normalize(q + vector));
            delta = endpoint-origin;
            h = quiver(origin(:,1),origin(:,2),delta(:,1),delta(:,2));
            h.LineWidth=2;
            h.Color = 0.1 * [1 1 1];
        end
        
    end
    
end