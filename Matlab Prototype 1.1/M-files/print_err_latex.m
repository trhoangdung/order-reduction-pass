function print_err_latex(sys_err,sys_bisim_err,fid)
% This function is used to print the experimental results in computing
% error bounds and time implement for different methods 
[me,ne] = size(sys_err);
err = sys_err.err;   % number of k (order)
[mO,nO] = size(err); % nO = number of system's outputs

for i = 1:ne
    j = i;
fprintf(fid,'& $%d$ \n',sys_err(i).order); % print order of output abstraction

fprintf(fid,'%%Prof Girard : bisimulation method \n'); % print result of bisimulation method
fprintf(fid,'& $%.2g$ & $%.2g$ \n',sys_bisim_err.err(i),sys_bisim_err.time(i));

fprintf(fid,'%%Zhi Han : simulation-based method \n'); % print result of Zhi Han's method
 if nO == 1 
   fprintf(fid,'&$%.2g$& $%.2g$ &{$5$} \n', sys_err(i).err(1).e_total_sim,sys_err(i).time_sim);  
 else
     str_delta1 = '';
     for k = 1:nO-1
       str_delta = sprintf('%.2g ',sys_err(i).err(k).e_total_sim);
       str_delta = [str_delta ' \\\\ '];
       str_delta1 = [str_delta1 str_delta];
     end
     str_delta2 = sprintf('%.2g',sys_err(i).err(nO).e_total_sim);
     str_delta = [str_delta1 str_delta2];
     
     str_time = sprintf('%.2g ',sys_err(i).time_sim);
     str = ['&$\\begin{pmatrix}' str_delta '\\end{pmatrix}$ & $' str_time '$ & {$5$} \n'];
     fprintf(fid,str);
 end


fprintf(fid,'%%Our method : mixed bound \n'); % print result of our optimization-based + simulation method

 if nO == 1
   fprintf(fid,'&$%.2g$ &$%.2g$ &$%.2g$ &$%.2g$ \n',sys_err(i).err(1).e1_opt,sys_err(i).err(1).e2_sim,sys_err(i).err(1).e_total_mix,sys_err(i).time_mix);  
 else
     
     str_delta1 = '';
     str_e10 = '';
     str_e20 = '';
     for k = 1:nO-1
       str_e1 = sprintf('%.2g ',sys_err(i).err(k).e1_opt);
       str_e1 = [str_e1 ' \\\\ '];
       str_e10 = [str_e10 str_e1];
       str_e2 = sprintf('%.2g ',sys_err(i).err(k).e2_sim);
       str_e2 = [str_e2 ' \\\\ '];
       str_e20 = [str_e20 str_e2];
       str_delta = sprintf('%.2g ',sys_err(i).err(k).e_total_mix);
       str_delta = [str_delta ' \\\\ '];
       str_delta1 = [str_delta1 str_delta];
       
     end
       str_e11 = sprintf('%.2g',sys_err(i).err(nO).e1_opt);
       str_e22 = sprintf('%.2g',sys_err(i).err(nO).e2_sim);
       str_delta2 = sprintf('%.2g',sys_err(i).err(nO).e_total_mix);
       
     str_e1 = [str_e10 str_e11];
     str_e2 = [str_e20 str_e22];
     str_delta = [str_delta1 str_delta2];
     
     str_e1 = ['&$\\begin{pmatrix}' str_e1 ' \\end{pmatrix}$ '];
     str_e2 = ['&$\\begin{pmatrix}' str_e2 ' \\end{pmatrix}$ '];
     str_delta = ['&$\\begin{pmatrix}' str_delta ' \\end{pmatrix}$ '];
     str_time  = sprintf('& $%.2g$ \n',sys_err(i).time_opt);
     
     str = [str_e1 str_e2 str_delta str_time];     
     fprintf(fid,str);
 end


fprintf(fid,'%%Our method : theoretical bound \n'); % print the theoretical results of our method

 if nO == 1
   fprintf(fid,'&$%.2g$ &$%.2g$ &$%.2g$ &$%.2g$ \\\\ \n',sys_err(i).err(1).e1_theo,sys_err(i).err(1).e2_theo,sys_err(i).err(1).e_total_theo,sys_err(i).time_theo);  
 else
     str_e10 = '';
     str_e20 = '';
     str_delta0 = '';
     for k = 1:nO-1
       str_e1 = sprintf('%.2g ',sys_err(i).err(k).e1_theo);
       str_e1 = [str_e1 ' \\\\ '];
       str_e10 = [str_e10 str_e1];
       str_e2 = sprintf('%.2g ',sys_err(i).err(k).e2_theo);
       str_e2 = [str_e2 ' \\\\ '];
       str_e20 = [str_e20 str_e2];
       str_delta = sprintf('%.2g ',sys_err(i).err(k).e_total_theo);
       str_delta = [str_delta ' \\\\ '];
       str_delta0 = [str_delta0 str_delta];
       
     end
       str_e11 = sprintf('%.2g',sys_err(i).err(nO).e1_theo);
       str_e22 = sprintf('%.2g',sys_err(i).err(nO).e2_theo);
       str_delta2 = sprintf('%.2g',sys_err(i).err(nO).e_total_theo);
       
     str_e1 = [str_e10 str_e11];
     str_e2 = [str_e20 str_e22];
     str_delta = [str_delta0 str_delta2];
     
     str_e1 = ['&$\\begin{pmatrix}' str_e1 ' \\end{pmatrix}$ '];
     str_e2 = ['&$\\begin{pmatrix}' str_e2 ' \\end{pmatrix}$ '];
     str_delta = ['&$\\begin{pmatrix}' str_delta ' \\end{pmatrix}$ '];
     str_time  = sprintf('& $%.2g$ ',sys_err(i).time_theo);
     str_time = [str_time ' \\\\ \n'];
     
     str = [str_e1 str_e2 str_delta str_time];     
     fprintf(fid,str);
 end
 
   if j == ne
      fprintf(fid,'\\hline \n\n'); 
   else 
      fprintf(fid,'\\cline{2-15} \n');
   end

end

end

