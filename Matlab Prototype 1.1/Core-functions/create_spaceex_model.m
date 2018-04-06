function create_spaceex_model(file_name,sys,stoptime)
% This function is use to create SpaceEx model for LTI system 
% this function generate xml file and cfg
% Author: Hoang-Dung Tran Feb-11-2016
Create_SpaceEx_Cfg_file(file_name,sys,stoptime);
Create_SpaceEx_xml_file(file_name,sys);

end



function Create_SpaceEx_Cfg_file(file_name,sys,stoptime)
% This function creates the configuration file for SpaceEx format
% author : Hoang-Dung Tran

x0_lb = sys.x0_lb;
x0_ub = sys.x0_ub;
u_lb = sys.u_lb;
u_ub = sys.u_ub;
[y0_lb,y0_ub] = Init_output_bound(sys);

[mX,nX] = size(x0_lb);
[mY,nY] = size(y0_lb); 
[mU,nU] = size(u_lb);
x_str = '';
for i=1:mX
x_spec = 'x%d >= %9.7f & x%d <= %9.7f & ';
x_str = [x_str sprintf(x_spec,i,x0_lb(i),i,x0_ub(i))]; % initial condition of state variables
end

y_str = '';
for i=1:mY
y_spec = 'y%d >= %9.7f & y%d <= %9.7f & ';
y_str = [y_str sprintf(y_spec,i,y0_lb(i),i,y0_ub(i))]; % initial condition of output
end

u_str = '';
for i=1:mU
u_spec = 'u%d >= %9.7f & u%d <= %9.7f & ';
%u_spec = 'u%d == %2.1f &';
u_str = [u_str sprintf(u_spec,i,u_lb(i),i,u_ub(i))]; % input constraint
end

t_str = sprintf('t==0 & stoptime == %4.2f"\n',stoptime); % initial condition of time variable
init_str = ['initially = " ' x_str y_str u_str t_str];

% create configuration file
file_name = [file_name '.cfg'];
fid = fopen(file_name,'w');
fprintf(fid,'# analysis option \n');
fprintf(fid,'system = "sys"\n');
fprintf(fid,init_str);
fprintf(fid,'scenario = "supp"\n');
fprintf(fid,'directions = "box"\n');
fprintf(fid,'sampling-time = 0.001\n');
fprintf(fid,'time-horizon = %4.2f\n',stoptime);
fprintf(fid,'iter-max = 10\n');

y_str = '';
for i=1:mY-1
y_spec = 'y%d, ';
y_str = [y_str sprintf(y_spec,i)]; % initial condition of output
end
y_str = [y_str sprintf('y%d',mY)];
y_str = ['output-variables = "t, ' y_str '"\n'];
fprintf(fid,y_str);
%fprintf(fid,'output-variables = "t,y"\n');
fprintf(fid,'output-format = "GEN"\n');
fprintf(fid,'rel-err = 1.0e-8\n');
fprintf(fid,'abs-err = 1.0e-12\n');

end

function Create_SpaceEx_xml_file(file_name,sys)
%This function is used to generate automatically the SpaceEx xml file for % Linear system
%The SpaceEx format is then used for reachability analysis
% author : Hoang-Dung

file_name = [file_name '.xml'];
fid = fopen(file_name,'w');
fprintf(fid,'<?xml version="1.0" encoding="iso-8859-1"?>\n');
fprintf(fid,'<sspaceex xmlns="http://www-verimag.imag.fr/xml-namespaces/sspaceex" version="0.2" math="SpaceEx">\n');
fprintf(fid,'  <component id="core_component">\n');

[mA,nA] = size(sys.sys.a);
[mB,nB] = size(sys.sys.b);
[mC,nC] = size(sys.sys.c);

% declare x 
x_paras = cell(mA,1); 
x_names = cell(mA,1);
types = sprintf('type="real" ');
locals = sprintf('local="false" ');
d1s = sprintf('d1="1" ');
d2s = sprintf('d2="1" ');
dynamics = sprintf('dynamics="any" ');

for i=1:mA
    x_names{i} = sprintf('name="x%d" ',i); 
    x_paras{i} = ['    <param ' x_names{i} types locals d1s d2s dynamics '/>\n'];
end
% declare y 

y_paras = cell(mC,1); 
y_names   = cell(mC,1);

for i = 1:mC
    y_names{i} = sprintf('name="y%d" ',i); 
    y_paras{i} = ['    <param ' y_names{i} types locals d1s d2s dynamics '/>\n'];
end 

% declare t
t_name = sprintf('name="t" ');
t_para = ['    <param ' t_name types locals d1s d2s dynamics '/>\n'];

% declare u 

u_paras = cell(nB,1);
u_names  = cell(nB,1); 
u_dynamic = sprintf('dynamics="const" ');

for i = 1:nB
    u_names{i} = sprintf('name="u%d" ',i); 
    u_paras{i} = ['    <param ' u_names{i} types locals d1s d2s u_dynamic '/>\n'];
end 

% declare stoptime 
stoptime_name = sprintf('name="stoptime" ');
stoptime_dynamic = sprintf('dynamics="const" ');
stoptime_para = ['    <param ' stoptime_name types locals d1s d2s stoptime_dynamic '/>\n'];

% prinf parameter 

for i = 1:mA 
    fprintf(fid,x_paras{i});
end

for i = 1:mC
    fprintf(fid,y_paras{i});
end 

fprintf(fid,t_para);

for i = 1:nB
     fprintf(fid,u_paras{i});
end

fprintf(fid,stoptime_para);
fprintf(fid,'    <location id="1" name="Model" x="362.0" y="430.0" width="426.0" height="610.0">\n');


% print invariant

[x_dot,y_invariant] = SpaceEx_dynamics_transform(sys);
x_dots = cell(mA,1);
y_invars = cell(mC,1); 

fprintf(fid,'     <invariant>t &lt;= stoptime \n &amp;');
for i = 1:mC-1
    y_invars{i} = char(y_invariant(i));
    y_invars_tempt = sprintf('y%d == ',i);
    y_invars{i} = [y_invars_tempt y_invars{i} '\n &amp;'];
    fprintf(fid,y_invars{i});
end 
    y_invars{mC} = char(y_invariant(mC)); % prinf final line of invariant
    y_invars_tempt = sprintf('y%d == ',mC); 
    y_invars{mC} = [y_invars_tempt y_invars{mC} '</invariant>\n'];
fprintf(fid,y_invars{mC});

% prinf flows 

fprintf(fid, '     <flow> ');

for i = 1:mA
    x_dots{i} = char(x_dot(i));
    x_dots_tempt = sprintf('x%d'' == ',i);
    x_dots{i} = [x_dots_tempt x_dots{i} '\n &amp;'];
    fprintf(fid,x_dots{i});
end

fprintf(fid, 't'' == 1</flow>\n'); 

% print location
fprintf(fid, '     </location>\n');
fprintf(fid,'   </component>');
fprintf(fid, '     <component id="sys">\n');

for i = 1:mA 
     x_paras{i} = ['    <param ' x_names{i} types locals d1s d2s dynamics 'controlled="true" />\n'];
     fprintf(fid,x_paras{i});
end

t_para = ['    <param ' t_name types locals d1s d2s dynamics 'controlled="true" />\n'];
stoptime_para = ['    <param ' stoptime_name types locals d1s d2s stoptime_dynamic 'controlled="true" />\n'];
fprintf(fid,t_para);
fprintf(fid,stoptime_para);

for i = 1:mC 
     y_paras{i} = ['    <param ' y_names{i} types locals d1s d2s dynamics 'controlled="true" />\n'];
     fprintf(fid,y_paras{i});
end

for i = 1:nB 
     u_paras{i} = ['    <param ' u_names{i} types locals d1s d2s u_dynamic 'controlled="true" />\n'];
     fprintf(fid,u_paras{i});
end

% print bind 

fprintf(fid,'     <bind component="core_component" as="model">\n');

map_x = cell(mA,1);
for i = 1:mA
    map_x{i} = sprintf('       <map key="x%d">x%d</map>\n',i,i);
    fprintf(fid,map_x{i});
end
map_t = sprintf('       <map key="t">t</map>\n');
map_stoptime = sprintf('       <map key="stoptime">stoptime</map>\n');

fprintf(fid,map_t);
fprintf(fid,map_stoptime); 

map_y = cell(mC,1);
for i = 1:mC
    map_y{i} = sprintf('       <map key="y%d">y%d</map>\n',i,i);
    fprintf(fid,map_y{i});
end

map_u = cell(nB,1);
for i = 1:nB
    map_u{i} = sprintf('       <map key="u%d">u%d</map>\n',i,i);
    fprintf(fid,map_u{i});
end

fprintf(fid,'    </bind>\n');
fprintf(fid,'  </component>\n');
fprintf(fid,'</sspaceex>');

end


function [ x_dot, y_dot ] = SpaceEx_dynamics_transform(sys)
% This function is used to calculate the dynamics of a linear system
% author: Hoang-Dung Tran

[mA,nA] = size(sys.sys.a);
[mB,nB] = size(sys.sys.b); 

x = sym('x', [mA 1]); % create state variable 
u = sym('u',[nB 1]); % create control input variable
x_dot = sys.sys.a*x + sys.sys.b*u; % calculate x_dot 
digits(5); % the number of significant digits is 5
x_dot = vpa(x_dot); % convert the coefficients to numeric constants
y_dot = sys.sys.c*x;  % *** this because we need to declare y as the invariant of the system in SpaceEx 
y_dot = vpa(y_dot);
end

function [y0_min, y0_max] = Init_output_bound(sys)
% This function is used to find the bound of initial output state
% It returns Y0 = [y0_min y0_max], we use thi bounds to generate SpacEx
% model of the system
  % 1) y0_min : the lower bound vector of the initial state of the output
  % 2) y0_max : the uper bound vector of the initial state of the output
  
 C = sys.sys.c; 
 lb = sys.x0_lb;
 ub = sys.x0_ub;
 [mC,nC] = size(C);
 y0_min = zeros(mC,1);
 y0_max = zeros(mC,1);
 
 for i = 1:mC
   C_temp = C(i,:);
   [x,y0_min(i)] = linprog(C_temp,[],[],[],[],lb,ub) ; 
   [x,y0_max(i)] = linprog(-C(i,:),[],[],[],[],lb,ub) ;
   y0_max(i) = -y0_max(i);      
 end
   Y0 = [y0_min y0_max];
end
