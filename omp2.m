function [hat_x,hat_y,Aug_t,pos_array,Psi_t,aug_y,T,Psi,wn] = omp2(w1,y,K,M,N,Phi)

%% OMP����ƥ��׷�ٷ��ع��ź�(��������L_1�������Ż�����)
%ƥ��׷�٣��ҵ�һ�����ǿ���ȥ���ռ�����������ص�С������������ȥ�������ǵ�����ӡ���������ظ�ֱ����������С����ǡ����͡��ռ������������ݡ�
Psi=fft(eye(N,N))/sqrt(N);                        %  ����Ҷ���任����
T=Phi*Psi;                                       %  �ָ�����(��������*�������任����)
hat_y=zeros(1,N);                                 %  ���ع�������(�任��)����   
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
Aug_t=[];                                         %  ��������(��ʼֵΪ�վ���)  
Psi_t=[];                                         %  ��������(��ʼֵΪ�վ���)
% Psi_t1=[];% Aug_t1=[]; 
r_n=y;                                            %  �в�ֵ
% A0=[];% w=[];
dn=y;
% wn1=zeros(1,1);
for times=1:K                                   %  ��������(�������������,�õ�������ΪK)
    for col=1:N                                 %  �ָ����������������
        product(col)=abs(T(:,col)'*r_n);          %  �ָ�������������Ͳв��ͶӰϵ��(�ڻ�ֵ) 
    end
    [val,pos]=max(product);                       %  ���ͶӰϵ����Ӧ��λ�ã����ҵ�һ�����ǿ���ȥ���ռ�����������ص�С��
    Aug_t=[Aug_t,T(:,pos)];                       %  ��������  
    Psi_t=[Psi_t,Psi(:,pos)];
%     Aug_t1=T(:,pos);        Psi_t1=Psi(:,pos);
%     A0=T(:,pos);
    T(:,pos)=zeros(M,1);                          %  ѡ�е������㣨ʵ����Ӧ��ȥ����Ϊ�˼��Ұ������㣩����������ȥ�������ǵ�����ӡ��
     pos_array(times)=pos;                         %  ��¼���ͶӰϵ����λ��
%     a1=Phi'*Phi;     a2=pinv(a1);     w1=a2*Phi';        
    aug_y=Psi_t'*w1*r_n;                         %  ���Ƽ��������С����
    wn1=LMS1(times,M,y,Aug_t);
%     wn1=[wn1',0]';
    r_n=r_n-Aug_t*aug_y;                         %  �в�
%     wn1=LMS1(times,M,r_n,Aug_t)
%   end
end
c=Aug_t'*Aug_t;
C=inv(c);
% C=cholesky(c,K); 
% xa=C*Aug_t'*y;                                    %��С���˷��ع�����ԭʼ�źŹ���ֵ
%% LMS
wn=LMS(wn1,K,M,dn,Aug_t);                            %LMS�����ع�����ԭʼ�źŹ���ֵ
% wn*2^9
%% �ָ�����
% wn2=[-4861,606-392i,607+392i,552+810i,551-807i,917+461i,918-461i,900-366i,900+367i,984+316i,985-318i,581-625i,582+624i,421+618i,421-620i,931-260i,931+261i,292+338i,292-338i,722-134i,721+136i,440+245i,440-244i,353-157i,352+157i,474+189i,475-188i,500+27i,500-26i,11+410i,11-411i,172+163i,172-163i,7+210i,8-210i,2]'/512;
% hat_y1(pos_array)=conj(wn2);
% hat_x1=real(Psi'*hat_y1.');
hat_y(pos_array)=conj(wn) ;                     %  �ع�����������
hat_x=real(Psi'*hat_y.') ;                     %  ���渵��Ҷ�任�ع��õ�ʱ���ź�

