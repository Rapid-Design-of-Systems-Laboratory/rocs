function [] = STT1_quad(n)

%clc; clear all; close all;

disp('Creating quad function')

%n=3; % number of points

x=sym('x',[1,n]);
x=sym(x,'real');
X=sym('X',[1,1]);
X=sym(X,'real');
y=sym('y',[1,n+1]);
y=sym(y,'real');

for j=1:n
for i=1:n
if i~=j
l(j,i)=(X-x(i))/(x(j)-x(i));
end
end
end

l=l+eye(n);

L=l(:,1);
for i=1:n
for j=1:n-1
L(i)=L(i)*l(i,j+1);
end
end

for i=1:n
w(i,1)=int(L(i),X);
end

W=w-subs(w,X,x(1));

quadsum=0;
for i=1:n
quadsum=quadsum+W(i)*y(i+1);
end

quadstart=y(1)-subs(quadsum,X,x(n));

quad=quadsum+quadstart;

% Create function and inputs
disp('Writing quad function')
name_numbers=num2str([n]);
function_name=['STT1_Quad',name_numbers];
fid=fopen([function_name,'.m'],'w');
fprintf(fid,['function Y = ',function_name,'(Data,X)\n\n']);
fprintf(fid,'n=%i;\n\n',n);

fprintf(fid,'[a, b]=size(X);\n');
fprintf(fid,'length=max([a b]);\n\n');

for i=1:n
   fprintf(fid,'x%i=Data(%i);\n',i,i); 
end
fprintf(fid,'\n');
for i=1:n+1
   fprintf(fid,'y%i=Data(%i);\n',i,i+n); 
end
fprintf(fid,'\n');

fprintf(fid,'for i=1:length\n');
fprintf(fid,'X1=X(i);\n');
text=char(quad);
fprintf(fid,'Y(i)=%s;\n',text);
fprintf(fid,'end\n\n');

fprintf(fid,'return\n\n');
fclose(fid);

return