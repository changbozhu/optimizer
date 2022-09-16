function closematlabpool 
  disp('Close MATLAB parallel Pool....');
   poolobj = gcp('nocreate');  
   delete(poolobj);  
end  