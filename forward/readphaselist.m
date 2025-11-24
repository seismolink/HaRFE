function [ray] = readphaselist(filename)
% function to read phase list for raysum
fid = fopen(filename,'rt');
n = 1;
m = 0;
while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
%     disp(tline)
    if mod(n,2) == 1
        m = m+1;
        ray(m).path = str2num(tline);
    else
        ray(m).type = str2num(tline);
    end
    n = n+1;
end
fclose(fid);
end