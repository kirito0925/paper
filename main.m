clc;
clear;%清理之前的工作区域数据
%zjc-version1.1
%待优化的角度：
%1.取出外部集合的p_chrom和m_chrom导入draw函数绘制甘特图
%2.引入集合覆盖率比较pareto解集的优劣
%3.计算切比雪夫聚合函数值时去除量纲的影响
%注意事项：HV值的计算是基于单一算法得到解进行坐标变换计算的,所以可能会出现单一算法的HV值比另一算法的HV小的情况,
%         此时我们需要将所有的算法的解集放置于同一个坐标系下（即将所有算法的解集合并起来后在坐标变换），才可以比较不同算法间HV值的差异
global  N H SH NM TM ps time;%变成全局变量避免传参频繁减少时间消耗
%N:工件数 
%H:每个工件对应的工序数 
%SH:总工序数 
%NM:每个工序可选的机器数
%TM:10
%time:工件i的第j道工序在机器t上加工的时间
global reference_point;
global max_point;
T=10;
ps=100;%种群大小
path='D:\myFFJSP\BFFJSP-main\BFFJSP-main\DATASET\DATASET\instance1\';
Data={'data1','data2','data3','data4','data5'};
xlsx_tail='.xlsx';
for file=1:4
   PATH=[path,Data{file},xlsx_tail];
[N,H,SH,NM,TM,time]=read_para(PATH,file);%读取数据
respath='result\';
tmp6='\';
respath=[respath,Data{file},tmp6];
%清除无关变量
clear tmp6

a=['mkdir ' respath];%创立写结果的文件夹
system(a);%利用Windows命令行命令执行dos命令
fprintf('%s %s\r\n','Calculating ',Data{file});
h_array=zeros(1,20);%用于存储20次独立运行的hv值
L_array=zeros(1,20);%用于存储每一次前沿中解集的个数
totalPF=[];%将20轮迭代的PF存入该变长数组，最后统计最大值和最小值，用于hv值的坐标系的统一
s={};%存放每轮pareto解集的工序染色体和机器染色体
for round=1:20%独立运行20次

[weights,neighbour]=init_weight(ps+1,2,T);    % 初始化权重向量矩阵 ps加1是为了保证权重向量不出现0   
[p_chrom,m_chrom]=initialNDS();%初始化种群随机生成工序码和机器码,非支配排序初始化
%由于三角模糊数是三维，所以必须创建cell存储
fitness=cell(ps,2);%存储最大完成时间和总的工作负载
for i=1:ps
    [fitness{i,1},fitness{i,2}]=fit(p_chrom(i,:),m_chrom(i,:));%fit函数求出每个种群的最大完工时间和总负载,存储的是模糊数的形式
end
obj=finalvalue(fitness);%将三角模糊数形式的最大完工时间和总负载基于模糊规则映射为函数值，生成种群大小的obj形式的矩阵
[reference_point,~]=min(obj);% 初始化参考点,参考点为所有个体的各个子目标函数值中的最小值组成，参考点就是每代种群finalvalue中最小的最大完工时间和最小的总负载
[max_point,~]=max(obj);%初始化最大参考点，计算切比雪夫函数值用
F=zeros(1,ps);%计数器用于统计最后20%的迭代中陷入局部最优的情况
 AS=pareto(fitness);%pareto(fit)，fit函数求出最大完工时间和总负载，AS中存储种群中非支配等级最高的种群的序号(总范围1到ps),即外部集合EP
for i=1:200 %最大迭代次数
     fprintf('%s %s %d %s %d\r\n',Data{file},'round',round,'iter',i);%写入txt文件
      for j=1:ps
        [new_p1,new_m1,new_p2,new_m2]=evolution(p_chrom,m_chrom,neighbour,j,T);%交叉变异生产新解
        f1=cell(1,2);
        f2=cell(1,2);
        [f1{1,1},f1{1,2}]=fit(new_p1,new_m1);%计算出子代染色体的最大完工时间和总负载，存储模糊数形式
        [f2{1,1},f2{1,2}]=fit(new_p2,new_m2);
        zz1=finalvalue(f1);%将模糊数形式数据映射为函数值
        zz2=finalvalue(f2);
        reference_point=min(reference_point,zz1);%比较参考点是否需要更新
        reference_point=min(reference_point,zz2);
        max_point=max(max_point,zz1);%比较最大值是否需要更新
        max_point=max(max_point,zz2);
        nei=neighbour(j,:);%取出第j组种群的领域
        for t=1:T
            k=nei(t);%取邻居 更新邻域个体
            tmp{1,1}=fitness{k,1};
            tmp{1,2}=fitness{k,2};
            zzold=finalvalue(tmp);%计算该领域的函数值
            [p_chrom(k,:),m_chrom(k,:),flag]=updates(weights(k,:),p_chrom(k,:),m_chrom(k,:),new_p1,new_m1,zzold,zz1);  % 更新邻域解，注意传的参数
            if flag==1 %flag=1时，更新了领域个体
                fitness{k,1}=f1{1,1};%更新领域个体的两个目标函数值
                fitness{k,2}=f1{1,2};
            elseif(i>100&&flag==0)%迭代到后期开始统计
                F(k)=F(k)+1;
            end
        end
        
        for t=1:T
            k=nei(t);
            tmp{1,1}=fitness{k,1};
            tmp{1,2}=fitness{k,2};
            zzold=finalvalue(tmp);
            [p_chrom(k,:),m_chrom(k,:),flag]=updates(weights(k,:),p_chrom(k,:),m_chrom(k,:),new_p2,new_m2,zzold,zz2);  % 更新领域解，注意传的参数
            if flag==1%更新适应值 减少解码的适应值计算
                fitness{k,1}=f2{1,1};
                fitness{k,2}=f2{1,2};
            elseif(i>100&&flag==0)%迭代到50%
                F(k)=F(k)+1;
            end
        end
        if(i>100&&F(j)>3&&(~ismember(j,AS)))%前沿个体不进行局部搜索 迭代后期进行局部搜索
            tmp{1,1}=fitness{j,1};
            tmp{1,2}=fitness{j,2};
            zzold=finalvalue(tmp);
            [p_chrom(j,:),m_chrom(j,:),flag]=VNS(p_chrom(j,:),m_chrom(j,:),zzold,weights(j,:));%进行变邻域搜索，比较函数值确定是否需要更新
            if flag==0
                F(j)=F(j)+1;%如果没有更新解
            else
                F(j)=0;
                [fitness{j,1},fitness{j,2}]=fit(p_chrom(j,:),m_chrom(j,:));
                fprintf('%s\r\n','VNS sucess');
            end
        end

      end
      if i>99
        AS=pareto(fitness);%到100代再更新外部集合EP，后续使用基于邻域和外部集合的选择策略
      end
end
[AS]=pareto(fitness);%存放非支配解的序号

L=length(AS);%非支配解个数
%记录最优解集的工序染色体和工件染色体,并绘制甘特图
%需要注意的是,p1和m1里面会有重复解出现
p1=zeros(ps,SH);
m1=zeros(ps,SH);
for i=1:L 
    p1(i,:)=p_chrom(AS(i),:);
    m1(i,:)=m_chrom(AS(i),:);
end 
s{1,round}=p1;%将每轮迭代出的工序和机器染色体存入，便于后续绘制甘特图使用
s{2,round}=m1;
%draw(p1(1,:),m1(1,:));%调用draw函数绘制甘特图

obj=finalvalue(fitness);%将最终解集中的模糊数映射为函数值

for i=1:L
    newobj(i,:)=obj(AS(i),:);%基于外部集合存放的序号读取最优解集的函数值  
end  

newobj=unique(newobj,'rows');%删除重复行
newobj=unique(newobj,'rows');%删除重复行
[L,~]=size(newobj);%读取到pareto解集的个数
L_array(round)=L;%记录每一次前沿解集的个数
for cc=1:L
    totalPF=[totalPF;newobj(cc,:)];%存放30轮的最终pareto解集
end
end
fmin   = min(min(totalPF,[],1),zeros(1,2));%min(a,[],1)：返回每一列的最小值，即比较每个目标函数上的最小值，规避掉0的情况
fmax   = max(totalPF,[],1)*1.1;%返回每一列的最大值，乘以1.1可以规避掉归一化后最大值为1的情况
current_index=1;
for round=1:20
    newobj=[];%取出每一次的PF前沿
    endindex=current_index+L_array(round)-1;
    for i=current_index:endindex
        newobj=[newobj;totalPF(i,:)];
    end
    current_index=current_index+L_array(round);
    h_array(round)=myHV(newobj,fmin,fmax);%计算每一轮pareto解集的HV值
    tmp5=newobj';%将求得的前沿解集转置

    tmp1='res';%根据字符串编写出结果存储的路径
    tmp2=num2str(round);
    tmp3='.txt';
    resPATH=[respath tmp1 tmp2 tmp3];%写入的文件路径
    fout=fopen(resPATH,'w');
    fprintf(fout,'%d\r\n',h_array(round));%写入hv值并换行
    fprintf(fout,'%5.2f %6.3f\r\n',tmp5);%由于matlab是按照列存储和读取矩阵，则写入文件的时候也是按照列填写，因此要转置存入
    fclose(fout);
end

sum=0;
for i=1:20
    sum=sum+h_array(i);
end
sum=sum/20;
best_hv=max(h_array);%读取最大的HV值
best_index=find(h_array==best_hv);%找到30轮中最大HV值出现的轮数
L=length(best_index);
if L>1 %如果不止一个结果则选择点个数最多的结果作为最优的
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

tmp1='res_mean.txt';%写入最后的hv的均值 和最好的hv值的前沿面的位置
resPATH=[respath tmp1];
fout=fopen(resPATH,'w');
fprintf(fout,'%s %s %s\r\n','average_hv','best_PF_index','best_hv');%写入hv值并换行
fprintf(fout,'%f %d %f\r\n',sum,t,best_hv);
fclose(fout);
fprintf('%s %s\r\n','Finish ',Data{file});
end