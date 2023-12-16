function varargout = Calcfield(imsize)
% Input:
% imsize: [Nx, Ny, Nslc] for 3D, or [Nx, Ny] for 2D

% Dengrong Jiang, Johns Hopkins BME, Jan 2017

dim = numel(imsize);

% 3D
if dim == 3
    % dc term
    dc = ones(imsize);
    % linear terms
    [x, y, z] = meshgrid(linspace(-1,1,imsize(2)), linspace(-1,1,imsize(1)),...
        linspace(-1,1,imsize(3)));
    % quadratic terms
    xy = x.*y;
    xz = x.*z;
    yz = y.*z;
    x2 = x.*x;
    y2 = y.*y;
    z2 = z.*z;
    % cubic terms
    x3 = x2.*x;
    y3 = y2.*y;
    z3 = z2.*z;
    xyz = x.*y.*z;
    x2y = x2.*y;
    x2z = x2.*z;
    xy2 = x.*y2;    
    xz2 = x.*z2;
    y2z = y2.*z;
    yz2 = y.*z2;
    % output
    varargout{1} = dc;
    varargout{2} = x;
    varargout{3} = y;
    varargout{4} = z;
    varargout{5} = xy;
    varargout{6} = xz;
    varargout{7} = yz;
    varargout{8} = x2;
    varargout{9} = y2;
    varargout{10} = z2;
    varargout{11} = xy2;
    varargout{12} = x2y;
    varargout{13} = x2z;
    varargout{14} = xz2;
    varargout{15} = y2z;
    varargout{16} = yz2;
    varargout{17} = xyz;
    varargout{18} = x3;
    varargout{19} = y3;
    varargout{20} = z3;
elseif dim == 2
    % dc term
    dc = ones(imsize);
    % linear terms
    [x, y] = meshgrid(linspace(-1,1,imsize(2)), linspace(-1,1,imsize(1)));
    % quadratic terms
    xy = x.*y;
    x2 = x.*x;
    y2 = y.*y;
    % cubic terms
    x3 = x2.*x;
    x2y = x2.*y;
    xy2 = x.*y2;
    y3 = y2.*y;
    % output
    varargout{1} = dc;
    varargout{2} = x;
    varargout{3} = y;
    varargout{4} = xy;
    varargout{5} = x2;
    varargout{6} = y2;
    varargout{7} = x3;
    varargout{8} = y3;
    varargout{9} = xy2;
    varargout{10} = x2y;
else
    fprintf('not supported dimension! for now only support 2D and 3D\n');
end


