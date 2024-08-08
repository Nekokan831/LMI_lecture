
%LMI lecture1
clear all;
% LMI条件の通し番号 
lmi_num = 0;

% ファジィ変数の最大最小
% Mさんの卒論
Omega_kappa= [0.0 9.5];

Omega_chi= [1.2699 1.0];

Omega_phi = [-2.5087 -0.7005];

fuzzy_counter = 0;

%model 
A_c = cell(8,1); % 8行1列のセル配列を初期化
for i = 1:2
    for j = 1:2
        for k = 1:2
            fuzzy_counter = fuzzy_counter + 1;
            A_c{fuzzy_counter} = [0 Omega_kappa(i) 0 0 0 0;
                               -Omega_kappa(i) 0 Omega_chi(j) 0 0 0;
                                0 0 0 Omega_phi(k) 0 0;
                                0 0 0 0 1 0;
                                0 0 0 0 0 0;
                                0 0 0 0 0 0];
        end
    end
end

% A行列の出力
% for index = 1:fuzzy_counter
%     disp(['A{' num2str(index) '} = ']);
%     disp(A{index});
% end

setlmis([]);%LMI宣言
P=lmivar(1,[6 1]);%リアプノフ関数 対角以外0にするなら-1？

% まず，Pは正定なのでその定義
% -lmi_numは，一つ目の条件が0より大きいという意味
lmi_num = lmi_num + 1;
lmiterm([-lmi_num 1 1 P], 1, 1)

% -lmi_numは，lmi_num番目の条件が0より大きいという意味
% A{i}までで -PA_i を表し，sでシンボリックになるので，-PA_i - A_i^TPになる
for i = 1:8
    lmi_num = lmi_num + 1;
    lmiterm([-lmi_num 1 1 P], -1, A_c{i},'s')
end


%LMI条件
% lmiterm([-1 1 1 P],1,1);
% lmiterm([2 1 1 P],1,A{1},'s');
% lmiterm([3 1 1 P],1,A{2},'s');

%LMIを解く
LMISYS=getlmis;
[tmin,xfeas] = feasp(LMISYS);%解けた値を格納
%                           %tmin<0なら解けてる
%                           %help feasp 参照
%解けたPの呼び出し
Po = dec2mat(LMISYS,xfeas,P)
Poeig = eig(Po) %固有値による解チェック
