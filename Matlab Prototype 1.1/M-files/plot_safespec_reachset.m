


function plot_safespec_reachset(sys,option1,option2)

% option1 is choose to draw safe or unsafe region or both
% option2 is choose to fill color in the regions or not : 
% if option2 = 'line' : nofill 
% if option2 = 'solid' : fill
% 


if strcmp(option1,'S') % plot safe region
    S_flag = 1;
else
    S_flag = 0;
end

if strcmp(option1,'U') % plot unsafe region
   U_flag = 1;
else
    U_flag = 0;
end


if strcmp(option1,'SU') % plot both safe and unsafe region
   SU_flag = 1;
else
   SU_flag = 0;
end


S = sys.S;
U = sys.U;
[m,n] = size(S);
nf = 0;

for i = 1:n
    Si = S(i); %  safe region
    Ui = U(i); %  unsafe region    
    
    if ispolytope(Si)
        A = Si.matrix_A; 
        A = A(:,any(A,1)); % remove zero column in A
        [mA,nA] = size(A); % nA is equal to number of outputs in the safety specification
        % what is the dim? 
        for j = 1:nA-1
            for k = j+1:nA
                dim = [j;k];
                nf = nf + 1; 
                fig = figure(nf);                
                if ~isempty(Ui.matrix_A)&&(U_flag||SU_flag)                   
                   plot_safespec_polytope(Ui,dim,option2);
                   hold on;
                end
                if ~isempty(Si.matrix_A)&&(S_flag||SU_flag)
                plot_safespec_polytope(Si,dim,option2);
                end
            end
        end
    
    end
    
    if isellipsoid(Si)
        P = Si.matrix_P;
        P = P(:,any(P,1)); % remove zero column in P
       [mP,nP] = size(P); % nP is equal to number of outputs in the safety specification
                % what is the dim? 
        for j = 1:nP-1
            for k = j+1:nP
                dim = [j;k];
                nf = nf+1;
                fig = figure(nf);                
                if ~isempty(Ui.matrix_P)&&(U_flag||SU_flag)                    
                   plot_safespec_polytope(Ui,dim);
                   hold on;
                end
                if ~isempty(Si.matrix_P)&&(S_flag||SU_flag)
                    plot_safespec_ellipsoid(Ui,dim); 
                end
            end
        end
  
    end
      
end


end



function  plot_safespec_polytope(safe_spec,dim,option)
% this function is used to plot the safety specification of the original
% system and its output abstraction on specific pair of outputs.
SP = safe_spec;

if strcmp(SP.type,'polytope')
    OK_flag1 = 1; 
else
    OK_flag1 = 0;
end

if strcmp(SP.name,'safe')
    OK_flag2 = 1; 
else
    OK_flag2 = 0;
end

if strcmp(option,'line')
    line_flag = 1;
else 
    line_flag = 0;
end

if strcmp(option,'solid')
    solid_flag = 1;
else     
    solid_flag = 0;
end
if ~line_flag && ~solid_flag
    disp('error in choosing option2 for drawing safe or unsafe region');
end

if OK_flag1
    P = polytope(SP.matrix_A,-SP.matrix_B);
    Q = projection(P,dim);
    V = extreme(Q);
    if OK_flag2
        if  line_flag
          plot_p2p(V,'b');  
        elseif solid_flag
          plot(Q,'b');   
        end
       
    else
        if  line_flag
          plot_p2p(V,'r');  
        elseif solid_flag
          plot(Q,'r');   
        end
    end
    
else
    disp('Error!!! safety specification is not a polytope');
end

end

function  plot_safespec_ellipsoid(safe_spec,dim)
% this function is used to plot 2D ellipsoid (x-x0)^TE(x-x0) = 1
SP = safe_spec;
P = SP.matrix_P;
R = SP.radius;
ny = dim(1);
nx = dim(2); 

if  strcmp(SP.type,'ellipsoid')&&~isempty(P)&&(R~=0)
    OK_flag1 = 1;
else 
    OK_flag1 = 0;
end

if OK_flag1
c = [SP.center(ny);SP.center(nx)];
E = zeros(2,2);
E(1,1) = P(ny,ny);
E(1,2) = P(ny,nx);
E(2,1) = P(nx,ny);
E(2,2) = P(nx,nx);
E = E/R^2;
else 
  disp('Error!!! safety specification is not correct, check matrix P and radius R');  
end 

if strcmp(SP.name,'safe')
  mpt_plotellip(E,c,'b');
else
  mpt_plotellip(E,c,'r');
end
    
end

function TF = ispolytope(safe_spec)
if strcmp(safe_spec.type,'polytope')
    TF = 1;
else 
    TF = 0; 
end
end 

function TF = isellipsoid(safe_spec)
if strcmp(safe_spec.type,'ellipsoid')
    TF = 1;
else 
    TF = 0; 
end
end 

function plot_p2p(P,color)
% This function is used to plot the safety specification with polygonal
% shape

[mP,nP] = size(P); 
for i = 1:mP-1
    plot([P(i,1),P(i+1,1)],[P(i,2),P(i+1,2)],'color',color);
    hold on;
end
    plot([P(mP,1),P(1,1)],[P(mP,2),P(1,2)],'color',color);
    hold on;
end

function varargout = plot_2d_vertices(fname,varargin)
% plot_2d_vertices    Plot a sequence of polygons defined by their
%                     vertices.
%
% plot_2d_vertices(fname)
% Plots the contents of the file "fname", which must be a two-column
% list of vertices separated by empty lines:
%    x11 y11
%    x12 y12
%    ...
%    x1n1 y1n1
%
%    x21 y21
%    x22 y22
%    ...
%    x2n2 y2n2
%
%    ...
%
% Each sequence of vertices defines a polygon. When an empty line is 
% encountered, a new polygon is started.
% 
% plot_2d_vertices(fname,...)
% Passes "..." as options to the patch command that draws the polygons.
% E.g., to plot in red color:
%    plot_2d_vertices(fname,'r')
%
% H=plot_2d_vertices(...)
% Returns a vector of handles to the patch objects.
% Allows the manipulation of the patches, e.g., to remove the outline:
%    set(H,'LineStyle','none')
%

if nargin==0
    error('Filename not specified.');
end
if nargin>1 
    color=varargin{1:end};
else
    color='b';
end

fid=fopen(fname,'rt');
if fid~=-1
    X=[];
    Y=[];
    H=[];
    while ~feof(fid)
        tline = fgetl(fid);
        if (~strcmp(tline,'') && ~feof(fid))
            d = sscanf(tline,'%g');
            X=[X;d(1)];
            Y=[Y;d(2)];
        else
            %h=patch(X,Y,color);
            %h=patch(X,Y,color,'EdgeColor','none');
            h=patch(X,Y,color,'EdgeColor',color,'LineWidth',4);
            X=[];
            Y=[];
            H=[H;h];
        end 
    end
    fclose(fid);
    if nargout>0
        varargout(1)={H};
    end
else
    error('Error: Could not open file.')
end

end

function V = get_vertices(safe_spec)

if strcmp(safe_spec.type,'polytope')
    A = safe_spec.matrix_A;
    B = safe_spec.matrix_B; 
    A = A(:,any(A,1)); % remove zero column in A
    P = polytope(A,-B);
    V = extreme(P);
else
    disp('safety specification is not a polytope');
    V = 0;
end

end