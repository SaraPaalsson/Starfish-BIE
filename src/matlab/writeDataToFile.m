function [] = writeDataToFile(dataW, varstr, filename)
%Writes data to file. Variable names specified in varstr. Variables saved
%in struct dataW
fileID = fopen(filename,'w');

nbrvar = size(varstr,2);

% fprintf(fileID, '%s',varstr{:});
varformat = []; sizeformat = [];
sizevec = [];
for i=1:nbrvar
    varformat = [varformat '%s '];
    sizevec = [sizevec size(dataW.(varstr{i}),2)];
    sizeformat = [sizeformat '%d '];
end
varformat = [varformat ' \n']; sizeformat = [sizeformat '\n'];

fprintf(fileID, varformat, varstr{:});
fprintf(fileID, sizeformat, sizevec); 

for i=1:nbrvar
   dataformat = [];
   tmp = dataW.(varstr{i}); 
   for j=1:sizevec(i)
       dataformat = [dataformat '%12.15f '];  
   end
   dataformat = [dataformat '\n'];
   fprintf(fileID,dataformat,tmp);
end

fclose(fileID);
end