
%LMI lecture1
clear all;
% LMI条件の通し番号 
lmi_num = 0;

%モーメント対角成分のX軸回り
raw = 1.155; % ρ[kg/m^3]
S = 0.42; % [m^2]
omega = 2; % [m]
I_xx = 0.1580; % [kg^2]
V_a = [9 15];
delta_zero = 0.2393; %[rad]
M_c = cell(2,1);
for i = 1:2
    M_c{i} = (1/2) * raw * S * omega * V_a(i)^2 * ( 0.1659 *  delta_zero + 0.2578);
end

% Aのファジィ変数の最大最小
% Mさんの卒論
Omega_kappa= [0.0 9.5];

Omega_chi= [8.2699 14.0];

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

B_c = cell(2,1);
for i = 1:2
    B_c{i} = [-1 0 0;
       0 0 0 ;
       0 0 0 ;
       0 0 0 ;
       0 M_c{i}/I_xx 0 ;
       0 0 1];
end

% A行列の出力
% for index = 1:fuzzy_counter
%     disp(['A{' num2str(index) '} = ']);
%     disp(A{index});
% end

%% X,M1 の定義
[sA_x,sA_y] = size(A_c{1}); % A行列のサイズ取得
[sB_x,sB_y] = size(B_c{1}); % B行列のサイズ取得

setlmis([]);%LMI宣言
[X,X_num,X_ch] = lmivar(1,[sA_x 1]) % Xの定義　大きさは(sA_x)×(sA_x)
[Xc,Xc_num,Xc_ch] = lmivar(1,[1 1]) % Xcの定義　大きさは１

% lmivarでは，返り値として
% [行列ごとの番号，行列の各要素の最終番号，各要素が行列の形として出てくる]
%
%

% まず，Pは正定なのでその定義
% -lmi_numは，一つ目の条件が0より大きいという意味
lmi_num = lmi_num + 1;
lmiterm([-lmi_num 1 1 P], 1, 1)
% lmiterm([-lmi_num 1 1 X], 1, A{i}','s')
% lmiterm([-lmi_num 1 1 M], -1, B','s')

% -lmi_numは，lmi_num番目の条件が0より大きいという意味
% A{i}までで -PA_i を表し，sでシンボリックになるので，-PA_i - A_i^TPになる
for i = 1:16
    for j = i:16
        lmi_num = lmi_num + 1;
        lmiterm([-lmi_num 1 1 P], -1, A_c{i},'s')
    end
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
