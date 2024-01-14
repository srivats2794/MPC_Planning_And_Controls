function [matrix,status] = make_rand_PSD_matrix(size,weight)
% This func makes a random "size" sized positive semi-definite matrix
% Tune its eigen values using the weight
   matrix= rand(size);

   matrix= matrix+matrix';
   matrix= matrix+weight*eye(size);

   if issymmetric(matrix)
       if all(eig(matrix)>=0)
           status ="Successfull";
       else
           status = "Not PSD";
       end
   else
       status = "Not Symmetric";
   end
end

