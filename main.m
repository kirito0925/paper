clc;
clear;%����֮ǰ�Ĺ�����������
%zjc-version1.1
%���Ż��ĽǶȣ�
%1.ȡ���ⲿ���ϵ�p_chrom��m_chrom����draw�������Ƹ���ͼ
%2.���뼯�ϸ����ʱȽ�pareto�⼯������
%3.�����б�ѩ��ۺϺ���ֵʱȥ�����ٵ�Ӱ��
%ע�����HVֵ�ļ����ǻ��ڵ�һ�㷨�õ����������任�����,���Կ��ܻ���ֵ�һ�㷨��HVֵ����һ�㷨��HVС�����,
%         ��ʱ������Ҫ�����е��㷨�Ľ⼯������ͬһ������ϵ�£����������㷨�Ľ⼯�ϲ�������������任�����ſ��ԱȽϲ�ͬ�㷨��HVֵ�Ĳ���
global  N H SH NM TM ps time;%���ȫ�ֱ������⴫��Ƶ������ʱ������
%N:������ 
%H:ÿ��������Ӧ�Ĺ����� 
%SH:�ܹ����� 
%NM:ÿ�������ѡ�Ļ�����
%TM:10
%time:����i�ĵ�j�������ڻ���t�ϼӹ���ʱ��
global reference_point;
global max_point;
T=10;
ps=100;%��Ⱥ��С
path='D:\myFFJSP\BFFJSP-main\BFFJSP-main\DATASET\DATASET\instance1\';
Data={'data1','data2','data3','data4','data5'};
xlsx_tail='.xlsx';
for file=1:4
   PATH=[path,Data{file},xlsx_tail];
[N,H,SH,NM,TM,time]=read_para(PATH,file);%��ȡ����
respath='result\';
tmp6='\';
respath=[respath,Data{file},tmp6];
%����޹ر���
clear tmp6

a=['mkdir ' respath];%����д������ļ���
system(a);%����Windows����������ִ��dos����
fprintf('%s %s\r\n','Calculating ',Data{file});
h_array=zeros(1,20);%���ڴ洢20�ζ������е�hvֵ
L_array=zeros(1,20);%���ڴ洢ÿһ��ǰ���н⼯�ĸ���
totalPF=[];%��20�ֵ�����PF����ñ䳤���飬���ͳ�����ֵ����Сֵ������hvֵ������ϵ��ͳһ
s={};%���ÿ��pareto�⼯�Ĺ���Ⱦɫ��ͻ���Ⱦɫ��
for round=1:20%��������20��

[weights,neighbour]=init_weight(ps+1,2,T);    % ��ʼ��Ȩ���������� ps��1��Ϊ�˱�֤Ȩ������������0   
[p_chrom,m_chrom]=initialNDS();%��ʼ����Ⱥ������ɹ�����ͻ�����,��֧�������ʼ��
%��������ģ��������ά�����Ա��봴��cell�洢
fitness=cell(ps,2);%�洢������ʱ����ܵĹ�������
for i=1:ps
    [fitness{i,1},fitness{i,2}]=fit(p_chrom(i,:),m_chrom(i,:));%fit�������ÿ����Ⱥ������깤ʱ����ܸ���,�洢����ģ��������ʽ
end
obj=finalvalue(fitness);%������ģ������ʽ������깤ʱ����ܸ��ػ���ģ������ӳ��Ϊ����ֵ��������Ⱥ��С��obj��ʽ�ľ���
[reference_point,~]=min(obj);% ��ʼ���ο���,�ο���Ϊ���и���ĸ�����Ŀ�꺯��ֵ�е���Сֵ��ɣ��ο������ÿ����Ⱥfinalvalue����С������깤ʱ�����С���ܸ���
[max_point,~]=max(obj);%��ʼ�����ο��㣬�����б�ѩ����ֵ��
F=zeros(1,ps);%����������ͳ�����20%�ĵ���������ֲ����ŵ����
 AS=pareto(fitness);%pareto(fit)��fit�����������깤ʱ����ܸ��أ�AS�д洢��Ⱥ�з�֧��ȼ���ߵ���Ⱥ�����(�ܷ�Χ1��ps),���ⲿ����EP
for i=1:200 %����������
     fprintf('%s %s %d %s %d\r\n',Data{file},'round',round,'iter',i);%д��txt�ļ�
      for j=1:ps
        [new_p1,new_m1,new_p2,new_m2]=evolution(p_chrom,m_chrom,neighbour,j,T);%������������½�
        f1=cell(1,2);
        f2=cell(1,2);
        [f1{1,1},f1{1,2}]=fit(new_p1,new_m1);%������Ӵ�Ⱦɫ�������깤ʱ����ܸ��أ��洢ģ������ʽ
        [f2{1,1},f2{1,2}]=fit(new_p2,new_m2);
        zz1=finalvalue(f1);%��ģ������ʽ����ӳ��Ϊ����ֵ
        zz2=finalvalue(f2);
        reference_point=min(reference_point,zz1);%�Ƚϲο����Ƿ���Ҫ����
        reference_point=min(reference_point,zz2);
        max_point=max(max_point,zz1);%�Ƚ����ֵ�Ƿ���Ҫ����
        max_point=max(max_point,zz2);
        nei=neighbour(j,:);%ȡ����j����Ⱥ������
        for t=1:T
            k=nei(t);%ȡ�ھ� �����������
            tmp{1,1}=fitness{k,1};
            tmp{1,2}=fitness{k,2};
            zzold=finalvalue(tmp);%���������ĺ���ֵ
            [p_chrom(k,:),m_chrom(k,:),flag]=updates(weights(k,:),p_chrom(k,:),m_chrom(k,:),new_p1,new_m1,zzold,zz1);  % ��������⣬ע�⴫�Ĳ���
            if flag==1 %flag=1ʱ���������������
                fitness{k,1}=f1{1,1};%����������������Ŀ�꺯��ֵ
                fitness{k,2}=f1{1,2};
            elseif(i>100&&flag==0)%���������ڿ�ʼͳ��
                F(k)=F(k)+1;
            end
        end
        
        for t=1:T
            k=nei(t);
            tmp{1,1}=fitness{k,1};
            tmp{1,2}=fitness{k,2};
            zzold=finalvalue(tmp);
            [p_chrom(k,:),m_chrom(k,:),flag]=updates(weights(k,:),p_chrom(k,:),m_chrom(k,:),new_p2,new_m2,zzold,zz2);  % ��������⣬ע�⴫�Ĳ���
            if flag==1%������Ӧֵ ���ٽ������Ӧֵ����
                fitness{k,1}=f2{1,1};
                fitness{k,2}=f2{1,2};
            elseif(i>100&&flag==0)%������50%
                F(k)=F(k)+1;
            end
        end
        if(i>100&&F(j)>3&&(~ismember(j,AS)))%ǰ�ظ��岻���оֲ����� �������ڽ��оֲ�����
            tmp{1,1}=fitness{j,1};
            tmp{1,2}=fitness{j,2};
            zzold=finalvalue(tmp);
            [p_chrom(j,:),m_chrom(j,:),flag]=VNS(p_chrom(j,:),m_chrom(j,:),zzold,weights(j,:));%���б������������ȽϺ���ֵȷ���Ƿ���Ҫ����
            if flag==0
                F(j)=F(j)+1;%���û�и��½�
            else
                F(j)=0;
                [fitness{j,1},fitness{j,2}]=fit(p_chrom(j,:),m_chrom(j,:));
                fprintf('%s\r\n','VNS sucess');
            end
        end

      end
      if i>99
        AS=pareto(fitness);%��100���ٸ����ⲿ����EP������ʹ�û���������ⲿ���ϵ�ѡ�����
      end
end
[AS]=pareto(fitness);%��ŷ�֧�������

L=length(AS);%��֧������
%��¼���Ž⼯�Ĺ���Ⱦɫ��͹���Ⱦɫ��,�����Ƹ���ͼ
%��Ҫע�����,p1��m1��������ظ������
p1=zeros(ps,SH);
m1=zeros(ps,SH);
for i=1:L 
    p1(i,:)=p_chrom(AS(i),:);
    m1(i,:)=m_chrom(AS(i),:);
end 
s{1,round}=p1;%��ÿ�ֵ������Ĺ���ͻ���Ⱦɫ����룬���ں������Ƹ���ͼʹ��
s{2,round}=m1;
%draw(p1(1,:),m1(1,:));%����draw�������Ƹ���ͼ

obj=finalvalue(fitness);%�����ս⼯�е�ģ����ӳ��Ϊ����ֵ

for i=1:L
    newobj(i,:)=obj(AS(i),:);%�����ⲿ���ϴ�ŵ���Ŷ�ȡ���Ž⼯�ĺ���ֵ  
end  

newobj=unique(newobj,'rows');%ɾ���ظ���
newobj=unique(newobj,'rows');%ɾ���ظ���
[L,~]=size(newobj);%��ȡ��pareto�⼯�ĸ���
L_array(round)=L;%��¼ÿһ��ǰ�ؽ⼯�ĸ���
for cc=1:L
    totalPF=[totalPF;newobj(cc,:)];%���30�ֵ�����pareto�⼯
end
end
fmin   = min(min(totalPF,[],1),zeros(1,2));%min(a,[],1)������ÿһ�е���Сֵ�����Ƚ�ÿ��Ŀ�꺯���ϵ���Сֵ����ܵ�0�����
fmax   = max(totalPF,[],1)*1.1;%����ÿһ�е����ֵ������1.1���Թ�ܵ���һ�������ֵΪ1�����
current_index=1;
for round=1:20
    newobj=[];%ȡ��ÿһ�ε�PFǰ��
    endindex=current_index+L_array(round)-1;
    for i=current_index:endindex
        newobj=[newobj;totalPF(i,:)];
    end
    current_index=current_index+L_array(round);
    h_array(round)=myHV(newobj,fmin,fmax);%����ÿһ��pareto�⼯��HVֵ
    tmp5=newobj';%����õ�ǰ�ؽ⼯ת��

    tmp1='res';%�����ַ�����д������洢��·��
    tmp2=num2str(round);
    tmp3='.txt';
    resPATH=[respath tmp1 tmp2 tmp3];%д����ļ�·��
    fout=fopen(resPATH,'w');
    fprintf(fout,'%d\r\n',h_array(round));%д��hvֵ������
    fprintf(fout,'%5.2f %6.3f\r\n',tmp5);%����matlab�ǰ����д洢�Ͷ�ȡ������д���ļ���ʱ��Ҳ�ǰ�������д�����Ҫת�ô���
    fclose(fout);
end

sum=0;
for i=1:20
    sum=sum+h_array(i);
end
sum=sum/20;
best_hv=max(h_array);%��ȡ����HVֵ
best_index=find(h_array==best_hv);%�ҵ�30�������HVֵ���ֵ�����
L=length(best_index);
if L>1 %�����ֹһ�������ѡ���������Ľ����Ϊ���ŵ�
    t=best_index(1);
    round=L_array(t);
    for i=1:L
        if round<L_array(best_index(i))
            round=L_array(best_index(i));
            t=best_index(i);
        end
    end
else
    t=best_index;
end

tmp1='res_mean.txt';%д������hv�ľ�ֵ ����õ�hvֵ��ǰ�����λ��
resPATH=[respath tmp1];
fout=fopen(resPATH,'w');
fprintf(fout,'%s %s %s\r\n','average_hv','best_PF_index','best_hv');%д��hvֵ������
fprintf(fout,'%f %d %f\r\n',sum,t,best_hv);
fclose(fout);
fprintf('%s %s\r\n','Finish ',Data{file});
end