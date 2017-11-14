#!/usr/bin/octave -qf

% ------------------------ START OF INPUT PARAMETER REGION ------------------- %

% Foil geometry
c = 1;                % Geometric chord length
s = 2;                % Span (along y-axis)
alpha = 0.087266;     % Angle of attack (in radians)
NACA = [0 0 1 5];     % NACA 4-digit designation as a row vector;

% Surface resolution parameters
Ni = 100;            % Number of interpolation points along the foil

% ------------------------- END OF INPUT PARAMETER REGION -------------------- %


% ---------------------------------- LICENCE  -------------------------------- %
%                                                                              %
%     Copyrighted 2011, 2012 by Håkon Strandenes, hakostra@stud.ntnu.no        %
%                                                                              % 
%     This program is free software: you can redistribute it and/or modify     %
%     it under the terms of the GNU General Public License as published by     %
%     the Free Software Foundation, either version 3 of the License, or        %
%     (at your option) any later version.                                      %
%                                                                              %
%     This program is distributed in the hope that it will be useful,          %
%     but WITHOUT ANY WARRANTY; without even the implied warranty of           %
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            %
%     GNU General Public License for more details.                             %
%                                                                              %
%     You should have received a copy of the GNU General Public License        %
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.    %
% ---------------------------------------------------------------------------- %


% Create a vector with x-coordinates, camber and thickness
beta=linspace(0,pi,Ni);
x = c*(0.5*(1-cos(beta)));
y_c = zeros(size(x));
y_t = zeros(size(x));
theta = zeros(size(x));


% Values of m, p and t
m = NACA(1)/100;
p = NACA(2)/10;
t = (NACA(3)*10 + NACA(4))/100;


% Calculate thickness
% The upper expression will give the airfoil a finite thickness at the trailing
% edge, witch might cause trouble. The lower expression is corrected to give 
% zero thickness at the trailing edge, but the foil is strictly speaking no
% longer a proper NACA airfoil.
%
% See http://turbmodels.larc.nasa.gov/naca4412sep_val.html
%     http://en.wikipedia.org/wiki/NACA_airfoil

%y_t = (t*c/0.2) * (0.2969.*(x/c).^0.5 - 0.1260.*(x/c) - 0.3516.*(x/c).^2 + 0.2843.*(x/c).^3 - 0.1015.*(x/c).^4);
y_t = (t*c/0.2) * (0.2969.*(x/c).^0.5 - 0.1260.*(x/c) - 0.3516.*(x/c).^2 + 0.2843.*(x/c).^3 - 0.1036.*(x/c).^4);


% Calculate camber
if (p > 0)
  % Calculate camber
  y_c = y_c + (m.*x/p^2) .* (2*p - x/c) .* (x < p*c);
  y_c = y_c + (m.*(c-x)/(1-p)^2) .* (1 + x/c - 2*p) .* (x >= p*c);


  % Calculate theta-value
  theta = theta + atan( (m/p^2) * (2*p - 2*x/c) ) .* (x < p*c);
  theta = theta + atan( (m/(1-p)^2) * (-2*x/c + 2*p) ) .* (x >= p*c);
end


% Calculate coordinates of upper surface
Xu = x - y_t.*sin(theta);
Yu = y_c + y_t.*cos(theta);


% Calculate coordinates of lower surface
Xl = x + y_t.*sin(theta);
Yl = y_c - y_t.*cos(theta);


% Rotate foil to specified angle of attack
upper = [cos(alpha), sin(alpha); -sin(alpha), cos(alpha)] * [Xu ; Yu];
lower = [cos(alpha), sin(alpha); -sin(alpha), cos(alpha)] * [Xl ; Yl];


% Merge upper and lower surface (NB: Assume that the trailing edge is sharp)
% (see comments w.r.t. thickness calculation above)
X = [ upper(1,:) lower(1,Ni-1:-1:2) ];
Y = [ upper(2,:) lower(2,Ni-1:-1:2) ];
N = length(X);


% Triangulate the end surface
tri = [1, 2, N];
for i=2:Ni-1
  tri = [tri; i, i+1, N-i+2];
end
for i=Ni+1:N-1
  tri = [tri; i, i+1, N-i+2];
end


% Make it 3D
X = [X X];
Y = [Y Y];
Z = [(s/2)*ones(1,N) -(s/2)*ones(1,N)];


% Triangulate the second end surface
tri = [tri; tri(:,2)+N, tri(:,1)+N, tri(:,3)+N];


% Triangulate the top and bottom
for i=1:N-1
  tri = [tri; i, N+i, i+1];
end
tri = [tri; N, 2*N, 1];
for i=N+1:(2*N-1)
  tri = [tri; i, i+1, i-N+1];
end
tri = [tri; 2*N, N+1, 1];


% Open file
fo = fopen('airfoil.stl', 'w');


% Write file
fprintf(fo, 'solid airfoil\n');

for i=1:length(tri)
  % Calculate normal vector
  AB = [X(tri(i,2)) - X(tri(i,1)), Y(tri(i,2)) - Y(tri(i,1)), Z(tri(i,2)) - Z(tri(i,1))];
  AC = [X(tri(i,3)) - X(tri(i,1)), Y(tri(i,3)) - Y(tri(i,1)), Z(tri(i,3)) - Z(tri(i,1))];
  n = cross(AB, AC) ./ norm(cross(AB, AC));
  
  % Write facet
  fprintf(fo, '  facet normal %e %e %e\n', n);
  fprintf(fo, '    outer loop\n');
  fprintf(fo, '      vertex %e %e %e\n', X(tri(i,1)), Y(tri(i,1)), Z(tri(i,1)));
  fprintf(fo, '      vertex %e %e %e\n', X(tri(i,2)), Y(tri(i,2)), Z(tri(i,2)));
  fprintf(fo, '      vertex %e %e %e\n', X(tri(i,3)), Y(tri(i,3)), Z(tri(i,3)));
  fprintf(fo, '    endloop\n');
  fprintf(fo, '  endfacet\n');
end

fprintf(fo, 'endsolid airfoil\n');



% Close file
fclose(fo);
