%LMI lecture2
clear all;

%model
a=0.1;
b=-0.1;
A{1}=[1 5+a*5;0 10-10*b];
A{2}=[1 5-5*a;0 10+10*b];
B{1}=[1+a;1+b];
B{2}=[1-a;1-b];

setlmis([]);%LMI宣言
X=lmivar(1,[2 1]);%リアプノフ関数の逆数
M{1}=lmivar(2,[1 2]);%FX
M{2}=lmivar(2,[1 2]);

%LMI条件
lmiterm([-1 1 1 X],1,1);

for i=[1 2]
lmiterm([-(1+i) 1 1 X],-1,A{i}','s');
lmiterm([-(1+i) 1 1 M{i}],B{i},1,'s');
end

lmiterm([-4 1 1 X],-1,A{1}','s');
lmiterm([-4 1 1 M{1}],B{2},1,'s');
lmiterm([-4 1 1 X],-1,A{2}','s');
lmiterm([-4 1 1 M{2}],B{1},1,'s');

%LMIを解く
LMISYS=getlmis;
[tmin,xfeas]=feasp(LMISYS)
%解けたものを呼び出し
Po=dec2mat(LMISYS,xfeas,X)
Poeig=eig(Po)%固有値による解チェック
Mo1=X*dec2mat(LMISYS,xfeas,M{1})
Mo2=X*dec2mat(LMISYS,xfeas,M{2})


