clear all;

isnp = 250;
ndof = 300*100*7;

nmodes = 60;
ndeim = 60;

% number of snapshot partitions
kmax = 4;

%----------------------------------

for k = 1:kmax

%---------------------

file = strcat('../../step1/snapshots/snapshot_u',num2str(k));
disp(file); fflush(stdout);
A = load(file);
X = 0;

ii = 0;
for j = 1:isnp
for i = 1:ndof
 ii = ii + 1; 
 X(i,j) = A(ii);
end
end


[U,S,V] = svd(X,'econ');


file = strcat('pod_modes',num2str(k));
qm = U(1:ndof,1:nmodes);
dlmwrite(file,qm,'delimiter','\t');


%------
% init

if (k==1)

 disp('init'); fflush(stdout);

 file = '../../step1/snapshots/init';
 q0 = load(file);
 qtil0 = qm'*q0;

 file = 'qtil0';
 dlmwrite(file,qtil0,'delimiter','\t');

end

%------
% orthonormality check

qmTqm = qm'*qm;
s1=0; s2=0; for i=1:nmodes for j=1:nmodes if(i==j) s1=s1+qmTqm(i,j); else s2=s2+qmTqm(i,j); end; end; end;
disp(k); fflush(stdout);
disp(s1/nmodes); fflush(stdout);
disp(s2/nmodes/(nmodes-1)); fflush(stdout); 

%------

% energy 

e = 0;
for i = 1:nmodes
 e = e + S(i,i);
end

et = 0;
for i = 1:isnp
 et = et + S(i,i);
end

disp("energy= "); disp(e/et); fflush(stdout);

%-----------------

if(1 == 2)

%----------------

%%%%%%%%%% DEIM
disp("DEIM");

file = strcat('../../step1/snapshots/snapshot_S',num2str(k));
disp(file); fflush(stdout);
A = load(file);
X = 0;

ii = 0;
for j = 1:isnp
for i = 1:ndof
 ii = ii + 1; 
 X(i,j) = A(ii);
end
end

[U,S,V] = svd(X,'econ');

UF = U(1:ndof,1:ndeim);

z1 = 0;
maxval = 0.0;
for i = 1:ndof
if(abs(UF(i,1))>maxval)
maxval = abs(UF(i,1));
z1 = i;
end
end

phi = 0.0;
P = 0.0;
r = 0.0;
pvector = 0;

phi(1:ndof,1) = UF(1:ndof,1);
P(1:ndof,1) = 0.0; P(z1,1)=1.0;
pvector(1) = z1;

%--------------

for l = 2:ndeim

 PTphi = (P')*phi;
 PTphil = (P')*UF(1:ndof,l);
 alpha = inv(PTphi)*PTphil;
 r = UF(1:ndof,l) - phi*alpha;

 zl = 0;
 maxval = 0.0;
 for i = 1:ndof
 if(abs(r(i))>maxval)
 maxval = abs(r(i));
 zl = i;
 end
 end

 phi(1:ndof,l) = UF(1:ndof,l);
 P(1:ndof,l) = 0.0; P(zl,l)=1.0;
 pvector(l) = zl;

end

%-------------

B_deim = (qm')*UF*inv((P')*UF);

file = strcat('B_deim',num2str(k));
dlmwrite(file,B_deim,'delimiter','\t');

file = strcat('pvector',num2str(k));
dlmwrite(file,pvector,'delimiter','\t');

disp('---------------'); fflush(stdout);

%-----------------------------------

%-----------------

endif
%% 1 == 2

%----------------


end

%-------------------------------
