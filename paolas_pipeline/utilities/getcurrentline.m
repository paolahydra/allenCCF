function line=getcurrentline
s=dbstack;
if length(s)>1 %called from a script, function, or live block: return the calling line
    line=s(2).line;
else
    % probably called from selection: will not work
    line = 1;
end
end