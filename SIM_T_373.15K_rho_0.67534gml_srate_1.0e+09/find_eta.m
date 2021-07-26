% File to find shear viscosity from NEMD Simulations
% 2021.05.20 F.Fleckenstein

clc
clear all
close all

% Parameters to determine fit
check_ss_1=1.5;             % checks for steady state in T,p,pyz
check_ss_2=0.01;           % checks for steady state in T,p,pyz
check_ss_cons=5;            % checks consecutive time

check=["v_T","v_p","v_pyz"];
eta=["v_eta"];
file="thermo.NVT.1.dat";
file2="in.Bulk_NVT_sllod";

% Start
input={"variableNequal";"variabletmpequal";"variablerho_inequal";"variables_rateequal";"variabledtequal"};
fid=fopen(file2);
d=textscan(fid,'%s','delimiter','\n','Whitespace','');
fid=fclose(fid);
d=string(d{:});
d=strrep(d,' ','');
for i=1:length(d)
    for j=1:length(input)
    get_input1 = strfind(d(i),input{j});
    if get_input1 ~= 0
        input_string{j} = erase(d(i),input{j});
        get_input2 = strfind(input_string{j},"#");
        if get_input2 ~= 0
            input_string{j} = extractBefore(input_string{j},"#");
        end
    end
    input_header(j)=erase(input{j},"variable");
    input_header(j)=erase(input_header(j),"equal");
    end
end



fid=fopen(file);
c=textscan(fid,'%s','delimiter','\n','Whitespace','');
fid=fclose(fid);
c=string(c{:});

for i=1:length(c)
    newStr=extractBefore(c(i),2);
    find_start(i,:)=newStr;  
end
start=sum(count(find_start,'#'))+1;
header=extractAfter(c(start-1),"# ");
header=strsplit(header);
c_cut=c(start:end);

for i=1:length(c_cut)
    thermo(i,:)=strsplit(c_cut(i));
end

thermo=str2double(thermo);

k = find(strcmp(header,check(1)));

for i=1:length(check)
    k(i) = find(strcmp(header, check(i)));
end

check_timestep="TimeStep";
timestep=find(strcmp(header,check_timestep));
timestep=thermo(:,timestep);

k(end+1)= find(strcmp(header,eta));

[B_end,o_max,start_ss]= return_ss (k,thermo,timestep,check_ss_1,check_ss_2,check_ss_cons);

ss_all=max(o_max)+max(start_ss);

if ss_all < length(thermo)-check_ss_cons   
    for i=1:length(k)
        subplot(2,2,i);
        plot(timestep(ss_all:end),thermo(ss_all:end,k(i)))
        title(header(k(i)))
        saveas(gcf,'plots.png')
    end
end

sol={'T','p','pyz','shear_eta'};

for i=1:length(sol)
    sol{2,i}=mean(thermo(ss_all:end,k(i)));
end

sol{1,end+1}='std_eta';
sol{2,end}=std(thermo(ss_all:end,k(end)));
sol{1,end+1}='steady_state_step';
sol{2,end}=ss_all*(timestep(2)-timestep(1));

%write data 
T1 = cell2table(input_string);
T1.Properties.VariableNames = cellstr(input_header);
T2 = cell2table(sol(2:end,:));
T2.Properties.VariableNames = sol(1,:);
T=[T1 T2];
fileID = fopen('result.txt','w');
writetable(T,'result.txt')
fclose(fileID);









function [B_end,o_max,start_ss]= return_ss (k,thermo,timestep,check_ss_1,check_ss_2,check_ss_cons)
for i=1:(length(k)-1)
n=1;
o=0;
    while n>0 && o<length(thermo)
    n=0;
    consec=1;
    
    vec=thermo(:,k(i));
    mean_ss(i)=mean(vec);
    std_ss(i)=std(vec);
    for j=1:length(thermo)
        if thermo(j,k(i))<mean_ss(i)-check_ss_1*std_ss(i) || thermo(j,k(i))>mean_ss(i)+check_ss_1*std_ss(i)
            vec(j)=0;
            mean_ss(i)=mean(vec(j+1:end));
            std_ss(i)=std(vec(j+1:end));
        end
    end
    
    valid=find(vec);
    for m=1:(numel(valid)-1)
        if valid(m+1) == (valid(m)+1)
            consec=consec+1;
               if consec > check_ss_cons
                    start_ss(1,i)=valid(m-check_ss_cons+1);
                    break
               end
        else 
            consec = 1;
        end
    end
    
    Y=thermo(start_ss(i)+o:end,k(i));
    x=timestep(start_ss(i)+o:end);
    X=[ones(size(x)) x];


    B(i,:) = X\Y;
    
    B_test=B(i,2).*(timestep(2)-timestep(1));


    if B_test>check_ss_2 || B_test<-check_ss_2
            n=1;
            o=o+1;
    end


    end
    B_end(i)=B_test;
    o_max(i)=o;
end
end

