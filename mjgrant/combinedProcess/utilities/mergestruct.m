% Copyright (c) 2011, Igor Kaufman
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

function Res = mergestruct(A,B)
%% Recursively merges fields and subfields of structures A and B to result structure Res
% Simple recursive algorithm merges fields and subfields of two structures
%   Example:
%   A.field1=1;
%   A.field2.subfield1=1;
%   A.field2.subfield2=2;
% 
%   B.field1=1;
%   B.field2.subfield1=10;
%   B.field2.subfield3=30;
%   B.field3.subfield1=1;
% 
%   C=mergestruct(A,B);
%

Res=[];
if nargin>0
    Res=A;
end;
if nargin==1 || isstruct(B)==0
    return;
end;    
  fnb=fieldnames(B);
  
  for i=1:length(fnb)
     s=char(fnb(i));
     oldfield=[];
     if (isfield(A,s))
         oldfield=getfield(A,s);
     end    
     newfield=getfield(B,s);
     if isempty(oldfield) || isstruct(newfield)==0
       Res=setfield(Res,s,newfield);     
     else
       Res=setfield(Res,s,mergestruct(oldfield, newfield));  
     end    
  end    

end


