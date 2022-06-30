function [output]=string2boolean(string)
   if strcmp(string,'false')
     output = false;
   else
     output = true;
   end
end
