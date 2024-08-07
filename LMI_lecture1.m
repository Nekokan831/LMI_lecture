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
% ここでは，行列のタイプは２：シンメトリック
% 行列１は２：２×２
% 行列１のタイプは0：すべて同じ要素の対角成分だけになる
% 行列１は２×３
% 行列１のタイプは1：フルランク（でシンメトリック）
X_5 = lmivar(2, [2 3]);
% こうすることで，
% x_1 & 0 & 0 & 0 //
% 0 & x_1 & 0 & 0 //
% 0 & 0 & x_2 & x_4 //
%0 & 0 & x_4 & x_3
%みたいな形になる
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





















