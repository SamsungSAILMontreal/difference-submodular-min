function [W] = read_bilmes(file);
% read speech data obtained from Jeff Bilmes's group
x = textread(file,'%s');
tab = [];
for i=12:4:length(x);
    tabloc = [ str2double(x(i)), str2double(x(i+1)), str2double(x(i+2))];
    if (tabloc(1)~=1) && (tabloc(2)~=2)
    tab = [ tab; tabloc-2  ];
    end
end
p = max(tab(:,1));
ngroups = max(tab(:,2))-p;

W = logical(sparse(tab(:,1),tab(:,2)-p,1,p,ngroups));


