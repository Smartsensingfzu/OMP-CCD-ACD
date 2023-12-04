function [pos_array,aug_y,hat_x,hat_y] = omp1(y,K,M,N,Phi)

%% OMP正交匹配追踪法重构信号(本质上是L_1范数最优化问题)
%匹配追踪：找到一个其标记看上去与收集到的数据相关的小波；在数据中去除这个标记的所有印迹；不断重复直到我们能用小波标记“解释”收集到的所有数据。

Psi=fft(eye(N,N))/sqrt(N);                        %  傅里叶正变换矩阵
T=Phi*Psi;                                       %  恢复矩阵(测量矩阵*正交反变换矩阵)
hat_y=zeros(1,N);                                 %  待重构的谱域(变换域)向量    
hat_y1=zeros(1,N); 
Aug_t=[];                                         %  增量矩阵(初始值为空矩阵)
Phi_t=[];                                         %  增量矩阵(初始值为空矩阵)
Psi_t=[];                                         %  增量矩阵(初始值为空矩阵)
r_n=y;                                            %  残差值
a1=[];
a2=[];
A0=[];
for times=1:K                                     %  迭代次数(有噪声的情况下,该迭代次数为K)
    for col=1:N                                   %  恢复矩阵的所有列向量
        product(col)=abs(T(:,col)'*r_n);          %  恢复矩阵的列向量和残差的投影系数(内积值) 
    end
    [val,pos]=max(product);                       %  最大投影系数对应的位置，即找到一个其标记看上去与收集到的数据相关的小波
    Aug_t=[Aug_t,T(:,pos)];                       %  矩阵扩充 
    A0=T(:,pos);
%     pos
    T(:,pos)=zeros(M,1);                          %  选中的列置零（实质上应该去掉，为了简单我把它置零），在数据中去除这个标记的所有印迹
      c=Aug_t'*Aug_t;
   C=inv(c);
%        C=cholesky(c,times);
    aug_y=C*Aug_t'*y ;          %  最小二乘,使残差最小
%     Aug_t'*Aug_t
%     inv(Aug_t'*Aug_t)
    r_n=y-Aug_t*aug_y;                            %  残差
    pos_array(times)=pos;                         %  纪录最大投影系数的位置
end
% Aug_t'*Aug_t;
% xi=inv(Aug_t'*Aug_t)*Aug_t'*y;
% hat_y1(pos_array)=conj(wn);
% hat_x1=real(Psi'*hat_y1.');
% aug_y*512
hat_y(pos_array)=conj(aug_y);                           %  重构的谱域向量
hat_x=real(Psi'*hat_y.');                       %  做逆傅里叶变换重构得到时域信号

