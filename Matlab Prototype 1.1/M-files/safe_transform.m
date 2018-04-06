 function [SP,UP] = safe_transform(sys,selected_bound) % this is the main function for safety specification transformation
       
       SP_ori = sys.S; 
       [mO,nO] = size(SP_ori);
       SP_ori_temp = SP_ori;
       err_temp = selected_bound;
       SP = safe_spec_class;
       UP = safe_spec_class;
       for i=1:nO
           if (ischar(SP_ori(i).name))
              % check type of safety specification : polytope or ellipsoid 
              if (strcmp(SP_ori(i).type, 'polytope'))||(strcmp(SP_ori(i).type, 'ellipsoid')) 
              OK_Flag1 = 1;
              else
              error('undefine the type of safety specification');
              OK_Flag1 = 0;
              end 
       
              % check name of safety specification : safe or unsafe     
              if (strcmp(SP_ori(i).name, 'safe'))||(strcmp(SP_ori(i).name, 'unsafe')) 
              OK_Flag2 = 1;
              else     
              error('undefine the name of safety specification');  
              OK_Flag2 = 0;
              end
              
              % transform safety specification
              if (strcmp(SP_ori(i).type, 'polytope'))&& (OK_Flag1) && (OK_Flag2)
              [SP(i),UP(i)] = safe_transform_polytope(SP_ori_temp(i),err_temp); 
              elseif (strcmp(SP_ori(i).type, 'ellipsoid'))&& (OK_Flag1) && (OK_Flag2)
              [SP(i),UP(i)] = safe_transform_ellipsoid(SP_ori_temp(i),err_temp);
              end
           end   
       end
       
      end 
       
       %These are local functions that is used for safety specification transformation
       
       % this function is based on the Lemma 2 in the paper for DEDS2016
         
       
       function [SP,UP] = safe_transform_polytope(SP_ori,selected_bound)
       
            A = SP_ori.matrix_A; % get matrix A for transformation
            B = SP_ori.matrix_B; % get matrix B for transformation
            [mA,nA] = size(A);
            [mB,nB] = size(B); 
            [me,ne] = size(selected_bound); 
            Delta = zeros(mA,1); % create the vector Delta
            B_safe = zeros(mB,1); % vector for new safety specification   SP
            B_unsafe = zeros(mB,1); % vector for new unsafe specification UP
            
            if (isempty(A)) % check empty matrix A
                error('invalid matrix. Matrix A is empty')
            elseif (isempty(B)) % check if matrix B is empty or has more than one column
                error('invalid matrix. Matrix B is empty');
            elseif (isempty(selected_bound)) || (ne > 1) % check the error bound vector is empty or contain more than 1 column
                error('invalid error bound')
            elseif (nA~=me)||(mA~=mB)
                error('Matrix A, B and error bound vector are not consistent')
            else 
                for i = 1:mA
                    for j = 1:nA 
                      Delta(i) = Delta(i)+abs(A(i,j))*selected_bound(j);  
                    end
                end
                
                for i = 1:mB                  
                    B_safe(i) = B(i)+Delta(i); % calculate the new matrix B_safe and B_unsafe
                    B_unsafe(i) = B(i) - Delta(i); 
                end
            end
            
            % return the transform safety and unsafe specification for the
            % abstraction
            if (strcmp(SP_ori.name, 'safe'))
               SP_ori.name = 'safe';
               SP_ori.matrix_B = B_safe;
               SP = SP_ori; 
               SP_ori.name = 'unsafe';
               SP_ori.matrix_B = B_unsafe; 
               UP = SP_ori; 
            elseif (strcmp(SP_ori.name, 'unsafe')) 
               SP_ori.name = 'safe';
               SP_ori.matrix_B = [];
               SP = SP_ori; 
               SP_ori.name = 'unsafe';
               SP_ori.matrix_B = B_unsafe; 
               UP = SP_ori; 
            end
  
       end
        
   % this function is based on the Lemma 3 in the paper for DEDS2016
       function [SP,UP] = safe_transform_ellipsoid(SP_ori,selected_bound)
            SP = safe_spec_class;
            UP = safe_spec_class;
            P = SP_ori.matrix_P; % get matrix P for transformation
           % [G,p] = chol(P); % check P is positive or not
            [U,S,V] = svd(P); % singular value decoposition for matrix P , note that P is symetric => U^T = V
            
            R = SP_ori.radius; % get the radius R of the safety specification
            
            [mP,nP] = size(P);
            [me,ne] = size(selected_bound);  
            Delta = zeros(mP,1);
            Delta_R = 0; % 

            
            if (~issymmetric(P)) % check empty matrix A
                error('invalid matrix. Matrix P is not symmetric ')
           % elseif (p>0) 
           %     error('invalid matrix. Matrix P is not positive');
            elseif (isempty(selected_bound)) || (ne > 1) % check the error bound vector is empty or contain more than 1 column
                error('invalid error bound')
            elseif (nP~=me)
                error('Matrix P and error bound vector are not consistent')
            else 
                for i = 1:mP
                    for j = 1:nP 
                      Delta(i) = Delta(i)+abs(U(i,j))*selected_bound(j);  
                    end
                end
                
                for i = 1:mP                    
                    Delta_R = Delta_R + S(i,i)*Delta(i)^2;
                end
                
                    Delta_R = sqrt(Delta_R); % get the final value of Delta_R
                    
            end
            
            % return the transform safety and unsafe specification for the
            % abstraction
            if (strcmp(SP_ori.name, 'safe'))
               SP = SP_ori;
               if(R>Delta_R)
               SP.radius = R - Delta_R;
               else
               SP.radius = 'error';
               end
               UP = SP_ori;
               UP.name = 'unsafe';
               UP.radius = R + Delta_R;  
            elseif (strcmp(SP_ori.name, 'unsafe')) 
               SP = SP_ori; 
               SP.name = 'safe';
               SP.matrix_P = [];
               SP.radius = [];
               
               UP = SP_ori;
               UP.radius = R + Delta_R; 
            end
  
       end     