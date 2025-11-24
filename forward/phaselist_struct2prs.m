function [phstr] = phaselist_struct2prs(filename)
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
phstr = '[''';
for i = 1:length(ray)
    for j = 1:length(ray(i).path)
        phstr = [phstr num2str(ray(i).path(j)-1)];
        if ray(i).type(j) == 1
            phstr = [phstr 'P'];
        elseif ray(i).type(j) == 2
            phstr = [phstr 'S'];
        elseif ray(i).type(j) == 3
            phstr = [phstr 'T'];
        elseif ray(i).type(j) == 4
            phstr = [phstr 'p'];
        elseif ray(i).type(j) == 5
            phstr = [phstr 's'];
        elseif ray(i).type(j) == 6
            phstr = [phstr 't'];
        end
    end
    if i~=length(ray)
        phstr = [phstr ''','''];
    end
end
phstr = [phstr ''']'];
end