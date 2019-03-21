function [] = enhance_polyfit(n,d)

disp('Creating poly function')

%n=4; % number of points
%d=1; % STT order, 0 means just actual

x=sym('x',[1,n]);
x=sym(x,'real');
a=sym('a',[1,n*(d+1)]);
a=sym(a,'real');
y=sym('y',[1,n*(d+1)]);
y=sym(y,'real');

for i=1:n
    for j=1:n*(d+1)
        A(i,j,1)=x(i)^(j-1);
    end
end

for j=1:d
    for i=1:n
        A(i,:,j+1)=diff(A(i,:,j),x(i));
    end
end

for i=1:d+1
    M(1+n*(i-1):n*i,:)=A(:,:,i);
end

M;
B=y';

% sol=linsolve(M,B);

% Create function and inputs
disp('Writing poly function')
name_numbers=num2str([n*10+d]);
function_name=['Enhance_Polyfit',name_numbers];
fid=fopen([function_name,'.m'],'w');
fprintf(fid,['function Y = ',function_name,'(Data,X)\n\n']);
fprintf(fid,'n=%i;\n',n);
fprintf(fid,'d=%i;\n\n',d);

fprintf(fid,'[a, b]=size(X);\n');
fprintf(fid,'length=max([a b]);\n\n');

for i=1:n
   fprintf(fid,'x%i=Data(%i);\n',i,i); 
end
fprintf(fid,'\n');
for i=1:n*(d+1)
   fprintf(fid,'y(%i)=Data(%i);\n',i,i+n);
   % fprintf(fid,'y%i=Data(%i);\n',i,i+n); 
end
fprintf(fid,'\n');

fprintf(fid,'M=[');
for i=1:n*(d+1)
    fprintf(fid,'[');
    for j=1:n*(d+1)
        text=char(M(i,j));
        fprintf(fid,'%s,',text);
    end
    fprintf(fid,'];\n');
end
fprintf(fid,'];\n\n');

% fprintf(fid,'transpose(y);\n\n');
% fprintf(fid,'rank(M)\n\n');
fprintf(fid,'A=inv(M)*transpose(y);\n\n');

% for i=1:n*(d+1)
%     text=char(sol(i));
%     fprintf(fid,'a(%i)=%s;\n',i,text);
% end
% fprintf(fid,'\n');

fprintf(fid,'Y=zeros(length,1);\n\n');
fprintf(fid,'for i=1:n*(d+1)\n');
fprintf(fid,'\tY(:)=Y(:)+A(i)*X(:).^(i-1);\n');
fprintf(fid,'end\n\n');

fprintf(fid,'return\n\n');
fclose(fid);

return