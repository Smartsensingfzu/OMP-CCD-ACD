function [hat_x,hat_y,Aug_t,pos_array,Psi_t,aug_y,T,Psi,wn] = omp2(w1,y,K,M,N,Phi)

%% OMP正交匹配追踪法重构信号(本质上是L_1范数最优化问题)
%匹配追踪：找到一个其标记看上去与收集到的数据相关的小波；在数据中去除这个标记的所有印迹；不断重复直到我们能用小波标记“解释”收集到的所有数据。
Psi=fft(eye(N,N))/sqrt(N);                        %  傅里叶正变换矩阵
T=Phi*Psi;                                       %  恢复矩阵(测量矩阵*正交反变换矩阵)
hat_y=zeros(1,N);                                 %  待重构的谱域(变换域)向量   
% u1=0;
% [T1,fy1]=mapminmax(T,-1,1);
% for p1=1:1024
%    for p2=(p1+1):1024
%        sum1=0;
%        sum2=0;
%      for p3=1:M
%          sum1=sum1+T1(p3,p1)*conj(T1(p3,p2));
%          sum2=sum2+sqrt((T1(p3,p1)^2)*(T1(p3,p2)^2));
%      end
%      u1=[u1,abs(sum1/sum2)];
%    end  
% end
% u2=max(u1) 
hat_y1=zeros(1,N);                    
Aug_t=[];                                         %  增量矩阵(初始值为空矩阵)  
Psi_t=[];                                         %  增量矩阵(初始值为空矩阵)
% Psi_t1=[];% Aug_t1=[]; 
r_n=y;                                            %  残差值
% A0=[];% w=[];
dn=y;
% wn1=zeros(1,1);
for times=1:K                                   %  迭代次数(有噪声的情况下,该迭代次数为K)
    for col=1:N                                 %  恢复矩阵的所有列向量
        product(col)=abs(T(:,col)'*r_n);          %  恢复矩阵的列向量和残差的投影系数(内积值) 
    end
    [val,pos]=max(product);                       %  最大投影系数对应的位置，即找到一个其标记看上去与收集到的数据相关的小波
    Aug_t=[Aug_t,T(:,pos)];                       %  矩阵扩充  
    Psi_t=[Psi_t,Psi(:,pos)];
%     Aug_t1=T(:,pos);        Psi_t1=Psi(:,pos);
%     A0=T(:,pos);
    T(:,pos)=zeros(M,1);                          %  选中的列置零（实质上应该去掉，为了简单我把它置零），在数据中去除这个标记的所有印迹
     pos_array(times)=pos;                         %  纪录最大投影系数的位置
%     a1=Phi'*Phi;     a2=pinv(a1);     w1=a2*Phi';        
    aug_y=Psi_t'*w1*r_n;                         %  近似计算替代最小二乘
    wn1=LMS1(times,M,y,Aug_t);
%     wn1=[wn1',0]';
    r_n=r_n-Aug_t*aug_y;                         %  残差
%     wn1=LMS1(times,M,r_n,Aug_t)
%   end
end
c=Aug_t'*Aug_t;
C=inv(c);
% C=cholesky(c,K); 
% xa=C*Aug_t'*y;                                    %最小二乘法重构出的原始信号估计值
%% LMS
wn=LMS(wn1,K,M,dn,Aug_t);                            %LMS方法重构出的原始信号估计值
% wn*2^9
%% 恢复序列
% wn2=[-4861,606-392i,607+392i,552+810i,551-807i,917+461i,918-461i,900-366i,900+367i,984+316i,985-318i,581-625i,582+624i,421+618i,421-620i,931-260i,931+261i,292+338i,292-338i,722-134i,721+136i,440+245i,440-244i,353-157i,352+157i,474+189i,475-188i,500+27i,500-26i,11+410i,11-411i,172+163i,172-163i,7+210i,8-210i,2]'/512;
% hat_y1(pos_array)=conj(wn2);
% hat_x1=real(Psi'*hat_y1.');
hat_y(pos_array)=conj(wn) ;                     %  重构的谱域向量
hat_x=real(Psi'*hat_y.') ;                     %  做逆傅里叶变换重构得到时域信号

