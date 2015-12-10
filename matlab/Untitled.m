min =1;
for i=1:1:length(test3)
   if(test3(i)~=0)
       if(min>abs(test3(i)))
           min = abs(test3(i));
           pos=i;
       end
   end
end