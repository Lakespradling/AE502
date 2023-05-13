clc
clear
close all
format long g
load GGM03C.txt
cnm=GGM03C(:,3);
snm=GGM03C(:,4);
n1=GGM03C(:,1);
m1=GGM03C(:,2);
cnmp=zeros(size(cnm));
snmp=zeros(size(snm));
CS=zeros(361);
for k=1:size(GGM03C,1)
    m=m1(k);
    n=n1(k);
    if m==0
        coeff=sqrt(2*n+1)*sqrt(factorial(n-m)/factorial(n+m));
        cnmp(k)=coeff*cnm(k);
        CS(n+1,m+1)=cnmp(k);
        snmp(k)=coeff*snm(k);
    else
        coeff=sqrt(2*(2*n+1))*sqrt(factorial(n-m)/factorial(n+m));
        cnmp(k)=coeff*cnm(k);
        CS(n+1,m+1)=cnmp(k);
        snmp(k)=coeff*snm(k);
        CS(m,n+1)=snmp(k);
    end
end
CS=CS(1:71,1:71);
fid = fopen('CS.txt','w');
fprintf(fid,'%s','[');
fprintf(fid,'\n');
for i=1:size(CS,1)
    temp=CS(i,:);    
    save CS.txt temp -ascii -double -append
end
fid = fopen('CS.txt','a');
fprintf(fid,'%s ','];');
fclose(fid);

