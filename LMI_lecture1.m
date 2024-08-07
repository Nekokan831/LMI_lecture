%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 困ったときには
% doc lmivar
% でドキュメント参照のこと

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lmivar(行列のタイプ, [行列が何×何か　どういうタイプか])
% ここでは，行列のタイプは１：シンメトリック
% 行列は２：２×２
% 行列のタイプは１：フルランク（でシンメトリック）
x_1 = lmivar(1, [2 1]);
% こうすることで，
% x_1 & x_3 //
% x_3 & x_2
%みたいな形になる
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lmivar(行列のタイプ, [行列が何×何か　どういうタイプか])
% ここでは，行列のタイプは１：シンメトリック
% 行列は２：２×２
% 行列のタイプは-1：対角成分だけになる
X_2 = lmivar(1, [2 -1]);
% こうすることで，
% x_1 & 0 //
% 0 & x_2
%みたいな形になる
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lmivar(行列のタイプ, [行列が何×何か　どういうタイプか])
% ここでは，行列のタイプは１：シンメトリック
% 行列は２：２×２
% 行列のタイプは0：すべて同じ要素の対角成分だけになる
X_3 = lmivar(1, [2 0]);
% こうすることで，
% x_1 & 0 //
% 0 & x_1
%みたいな形になる
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lmivar(行列のタイプ, [行列が何×何か　どういうタイプか])
% [O O ; O O] とすることもできる．
%　これは，二つの行列を組み合わせたような形になる．
% ここでは，行列のタイプは１：シンメトリック
% 行列１は２：２×２
% 行列１のタイプは0：すべて同じ要素の対角成分だけになる
% 行列１は２：２×２
% 行列１のタイプは1：フルランク（でシンメトリック）
X_4 = lmivar(1, [2 0; 2 1]);
% こうすることで，
% x_1 & 0 & 0 & 0 //
% 0 & x_1 & 0 & 0 //
% 0 & 0 & x_2 & x_4 //
%0 & 0 & x_4 & x_3
%みたいな形になる
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lmivar(行列のタイプ, [行列が何×何か　どういうタイプか])
% [O O ; O O] とすることもできる．
%　これは，二つの行列を組み合わせたような形になる．
% ここでは，行列のタイプは２：長方形
% 行列１は２：２×３
X_5 = lmivar(2, [2 3]);
% こうすることで，
% x_1 & x_2 & x_3 //
% x_4 & x_5 & x_6
%みたいな形になる
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%LMI lecture1
clear all;

%model 
k=3.9; %parameter(under3.9)
A{1}=[0 1;-2 -1];
A{2}=[0 1;-2-k -1];

setlmis([]);%LMI宣言
P=lmivar(1,[2 1]);%リアプノフ関数

%LMI条件
lmiterm([-1 1 1 P],1,1);
lmiterm([2 1 1 P],1,A{1},'s');
lmiterm([3 1 1 P],1,A{2},'s');

%LMIを解く
LMISYS=getlmis;
[tmin,xfeas]=feasp(LMISYS)%解けた値を格納
                          %tmin<0なら解けてる
                          %help feasp 参照
%解けたPの呼び出し
Po=dec2mat(LMISYS,xfeas,P)
Poeig=eig(Po)%固有値による解チェック




















