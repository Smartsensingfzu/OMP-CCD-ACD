function [pos_array,aug_y,hat_x,hat_y] = omp1(y,K,M,N,Phi)

%% OMP����ƥ��׷�ٷ��ع��ź�(��������L_1�������Ż�����)
%ƥ��׷�٣��ҵ�һ�����ǿ���ȥ���ռ�����������ص�С������������ȥ�������ǵ�����ӡ���������ظ�ֱ����������С����ǡ����͡��ռ������������ݡ�

Psi=fft(eye(N,N))/sqrt(N);                        %  ����Ҷ���任����
T=Phi*Psi;                                       %  �ָ�����(��������*�������任����)
hat_y=zeros(1,N);                                 %  ���ع�������(�任��)����    
hat_y1=zeros(1,N); 
Aug_t=[];                                         %  ��������(��ʼֵΪ�վ���)
Phi_t=[];                                         %  ��������(��ʼֵΪ�վ���)
Psi_t=[];                                         %  ��������(��ʼֵΪ�վ���)
r_n=y;                                            %  �в�ֵ
a1=[];
a2=[];
A0=[];
for times=1:K                                     %  ��������(�������������,�õ�������ΪK)
    for col=1:N                                   %  �ָ����������������
        product(col)=abs(T(:,col)'*r_n);          %  �ָ�������������Ͳв��ͶӰϵ��(�ڻ�ֵ) 
    end
    [val,pos]=max(product);                       %  ���ͶӰϵ����Ӧ��λ�ã����ҵ�һ�����ǿ���ȥ���ռ�����������ص�С��
    Aug_t=[Aug_t,T(:,pos)];                       %  �������� 
    A0=T(:,pos);
%     pos
    T(:,pos)=zeros(M,1);                          %  ѡ�е������㣨ʵ����Ӧ��ȥ����Ϊ�˼��Ұ������㣩����������ȥ�������ǵ�����ӡ��
      c=Aug_t'*Aug_t;
   C=inv(c);
%        C=cholesky(c,times);
    aug_y=C*Aug_t'*y ;          %  ��С����,ʹ�в���С
%     Aug_t'*Aug_t
%     inv(Aug_t'*Aug_t)
    r_n=y-Aug_t*aug_y;                            %  �в�
    pos_array(times)=pos;                         %  ��¼���ͶӰϵ����λ��
end
% Aug_t'*Aug_t;
% xi=inv(Aug_t'*Aug_t)*Aug_t'*y;
% hat_y1(pos_array)=conj(wn);
% hat_x1=real(Psi'*hat_y1.');
% aug_y*512
hat_y(pos_array)=conj(aug_y);                           %  �ع�����������
hat_x=real(Psi'*hat_y.');                       %  ���渵��Ҷ�任�ع��õ�ʱ���ź�

