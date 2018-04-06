
% This mfile is used to print a latex table of the error bound and
% computing time 
% input 1: E class
% input 2: file name want to be saved with 

% Print table of experimental result
file_name = ['table_result' '.tex'];
fid = fopen(file_name,'w');
fprintf(fid,'\\begin{table*}[t] \n');
fprintf(fid,'\\centering \n');
fprintf(fid,'\\scriptsize\n');
%fprintf(fid,'\\tiny \n');
fprintf(fid,'\\tabcolsep=0.11cm \n');
fprintf(fid,'\\begin{tabular}{|p{1.5cm}|c|c|c|c|c|c|c|c|c|c|c|c|c|c|} \n');
fprintf(fid,'\\hline \n');
fprintf(fid,'\\multirow{2}{*}{\\textbf{Benchmark}} &\\multirow{2}{*}{\\textbf{k}} & \\multicolumn{2}{c|}{\\textbf{MATISSE~\\cite{girard2007approximate}}} & \\multicolumn{3}{c|}{\\textbf{Zhi Han~\\cite{han2004reachability}}} & \\multicolumn{4}{c|}{\\textbf{Mixed bound}} & \\multicolumn{4}{c|}{\\textbf{Theoretical bound}} \\\\ \n');
fprintf(fid,'\\cline{3-15} \n');
fprintf(fid,'&  & {$\\delta$} & $t$ & $\\delta$ & $t$ & $N$ & $ e_1 $ & $ e_2 $ & $\\delta$ &$t$ & $e_1 $ & $e_2$ & $\\delta$ &$t$ \\\\ \n');
fprintf(fid,'\\hline \n');

% print data 

% print motor control system data : 8 dimensions, 2 inputs, 2 outputs
[me,ne] = size(mcs_err);
fprintf(fid,'\\multirow{%d}{*}{\\parbox{1.5cm}{Motor control system}} \n',ne);
print_err_latex(mcs_err,mcs_bisim_err,fid);

% print Helicopter data : 28 dimensions, 6 inputs, 2 outputs
[me,ne] = size(heli_err);
fprintf(fid,'\\multirow{%d}{*}{\\parbox{1.5cm}{Helicopter}} \n',ne);
print_err_latex(heli_err,heli_bisim_err,fid);

% print Building model data : 48 dimensions, 1 inputs, 1 outputs
[me,ne] = size(build_err);
fprintf(fid,'\\multirow{%d}{*}{\\parbox{1.5cm}{Building model}} \n',ne);
print_err_latex(build_err,build_bisim_err,fid);

% print Partial different equation Pde data : 84 dimensions, 1 inputs, 1 outputs
[me,ne] = size(pde_err);
fprintf(fid,'\\multirow{%d}{*}{\\parbox{1.5cm}{Partial differential equation}} \n',ne);
print_err_latex(pde_err,pde_bisim_err,fid);

% print Internation space station ISS data : 270 dimensions, 3 inputs, 3 outputs
[me,ne] = size(iss_err);
fprintf(fid,'\\multirow{%d}{*}{\\parbox{1.5cm}{International space station}} \n',ne);
print_err_latex(iss_err,iss_bisim_err,fid);

% print FOM model data : 1006 dimensions, 1 inputs, 1 outputs
[me,ne] = size(fom_err);
fprintf(fid,'\\multirow{%d}{*}{\\parbox{1.5cm}{FOM model}} \n',ne);
print_err_latex(fom_err,fom_bisim_err,fid);

fprintf(fid,'\\end{tabular} \n');
fprintf(fid,' \\caption{The error bounds and computing times obtained from different methods on different benchmarks in which: $k$ is the dimension of the output abstraction, $\\delta$ is total error bound, $e_1$ is the zero input response error, $e_2$ is the zero state response error, $t$ is the error computing time (in second) and $N$ is the number of simulations.} \n');
fprintf(fid,'\\tablabel{experiment_test} \n');
fprintf(fid,'\\end{table*} \n');




