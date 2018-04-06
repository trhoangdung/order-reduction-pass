function error_print(Err,file_name)

% This function is used to print the result of the error to a .txt file
% input 1: Error class
% input 2: file name want to be saved with 

[m,n] = size(Err);

file_name = [file_name '.txt'];
fid = fopen(file_name,'w');
fprintf(fid,'-----------The information of the errors between original system and reduced system------------ \n\n');

err_string = ''; 
for i = 1:n
    fprintf(fid,'--------------ERROR %d-------------- \n\n',i);
    
    fprintf(fid,'e1 theory bound: '); 
    err_string = num2str(Err(i).e1_theo);
    err_string = [err_string '\n\n'];
    fprintf(fid,err_string);
    
    fprintf(fid,'e1 simulation bound: ')
    err_string = num2str(Err(i).e1_sim);
    err_string = [err_string '\n\n'];
    fprintf(fid,err_string);
    
    fprintf(fid,'e1 optimized bound: ')
    err_string = num2str(Err(i).e1_opt);
    err_string = [err_string '\n\n'];
    fprintf(fid,err_string);
    
    fprintf(fid,'e2 theo bound: ')
    err_string = num2str(Err(i).e2_theo);
    err_string = [err_string '\n\n'];
    fprintf(fid,err_string);
    
    fprintf(fid,'e2 simulation bound: ')
    err_string = num2str(Err(i).e2_sim);
    err_string = [err_string '\n\n'];
    fprintf(fid,err_string);
    
    fprintf(fid,'e2 optimized bound: ')
    err_string = num2str(Err(i).e2_opt);
    err_string = [err_string '\n\n'];
    fprintf(fid,err_string);
 
    fprintf(fid,'total error theory bound: ')
    err_string = num2str(Err(i).e_total_theo);
    err_string = [err_string '\n\n'];
    fprintf(fid,err_string);
    
    fprintf(fid,'total error simulation bound: ')
    err_string = num2str(Err(i).e_total_sim);
    err_string = [err_string '\n\n'];
    fprintf(fid,err_string);
    
    fprintf(fid,'total error mixed bound: ')
    err_string = num2str(Err(i).e_total_mix);
    err_string = [err_string '\n\n'];
    fprintf(fid,err_string);
    
    
    fprintf(fid,'computing time of mixed bound: ')
    err_string = num2str(Err(i).time_opt);
    err_string = [err_string '\n\n'];
    fprintf(fid,err_string);
    
        
    fprintf(fid,'computing time of theoretical bound: ')
    err_string = num2str(Err(i).time_theo);
    err_string = [err_string '\n\n'];
    fprintf(fid,err_string);
       
end

end

