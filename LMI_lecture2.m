
%LMI lecture1
clear all;

% ファジィ変数の最大最小
% Mさんの卒論
Omega_kappa= [0.0 9.5];

Omega_chi= [8.2699 14.0];

Omega_phi = [-2.5087 -0.7005];

fuzzy_counter = 0;

%model 
for i = 1:2
    for j = 1:2
        for k = 1:2
            fuzzy_counter += 1;
            A{fuzzy_counter}=[0 Omega_kappa[i] 0 0 0;
                            -Omega_kappa[i] 0 Omega_chi[j] 0 0;
                            0 0 0 Omega_phi[k] 0;
                            0 0 0 0 1;
                            0 0 0 0 0];
        end
    end
end

% setlmis([]);%LMI宣言
% P=lmivar(1,[2 1]);%リアプノフ関数

% %LMI条件
% lmiterm([-1 1 1 P],1,1);
% lmiterm([2 1 1 P],1,A{1},'s');
% lmiterm([3 1 1 P],1,A{2},'s');

% %LMIを解く
% LMISYS=getlmis;
% [tmin,xfeas]=feasp(LMISYS)%解けた値を格納
%                           %tmin<0なら解けてる
%                           %help feasp 参照
% %解けたPの呼び出し
% Po=dec2mat(LMISYS,xfeas,P)
% Poeig=eig(Po)%固有値による解チェック
