function F = newton_euler_lean(X)

global m_foot_L m_leg_L m_thigh_L m_pelvis_L m_arm_L m_torso ...
    m_foot_R m_leg_R m_thigh_R m_pelvis_R m_arm_R m_pelvis m_head_neck ...
    m_torso_arms g_x g_y g_z m_total head_scale num_rows ... 
    foot_L_cg_x foot_L_cg_y foot_L_cg_z leg_L_cg_x leg_L_cg_y leg_L_cg_z ...
    thigh_L_cg_x thigh_L_cg_y thigh_L_cg_z arm_L_cg_x arm_L_cg_y arm_L_cg_z...
    torso_L_cg_x torso_L_cg_y torso_L_cg_z pelvis_L_cg_x pelvis_L_cg_y ...
    pelvis_L_cg_z torso_cg_x torso_cg_y torso_cg_z torso_R_cg_x torso_R_cg_y...
    torso_R_cg_z head_neck_cg_x head_neck_cg_y head_neck_cg_z foot_R_cg_x...
    foot_R_cg_y foot_R_cg_z leg_R_cg_x leg_R_cg_y leg_R_cg_z thigh_R_cg_x ...
    thigh_R_cg_y thigh_R_cg_z arm_R_cg_x arm_R_cg_y arm_R_cg_z pelvis_R_cg_x...
    pelvis_R_cg_y pelvis_R_cg_z   r_1_L_y r_1_L_z r_1_R_y r_1_R_z r_2_L_x r_2_L_y r_2_L_z r_2_R_x r_2_R_y r_2_R_z...
    r_3_L_x r_3_L_y r_3_L_z r_3_R_x r_3_R_y r_3_R_z r_4_L_x r_4_L_y r_4_L_z...
    r_4_R_x r_4_R_y r_4_R_z r_5_L_x r_5_L_y r_5_L_z r_5_R_x r_5_R_y r_5_R_z...
    r_6_L_x r_6_L_y r_6_L_z r_6_R_x r_6_R_y r_6_R_z r_7_L_x r_7_L_y r_7_L_z ...
    r_7_R_x r_7_R_y r_7_R_z r_8_x r_8_y r_8_z r_9_x r_9_y r_9_z r_10_x r_10_y r_10_z
%Equations of motion inside function
%for loop after all constants have been defined in run file.
%r_3_L etc to be defined in run file for loop

F_GRF_L_x = X(1);
F_GRF_L_y = X(2);
F_GRF_L_z = X(3);
F_GRF_R_x = X(4);
F_GRF_R_y = X(5);
F_GRF_R_z = X(6);
F_ANKLE_L_x = X(7);
F_ANKLE_L_y = X(8);
F_ANKLE_L_z = X(9);
F_ANKLE_R_x = X(10);
F_ANKLE_R_y = X(11);
F_ANKLE_R_z = X(12);
F_KNEE_L_x = X(13);
F_KNEE_L_y = X(14);
F_KNEE_L_z = X(15);
F_KNEE_R_x = X(16);
F_KNEE_R_y = X(17);
F_KNEE_R_z = X(18);
F_HIP_L_x = X(19);
F_HIP_L_y = X(20);
F_HIP_L_z = X(21);
F_HIP_R_x = X(22);
F_HIP_R_y = X(23);
F_HIP_R_z = X(24);
F_PELVIS_x = X(25);
F_PELVIS_y = X(26);
F_PELVIS_z = X(27);
F_HEAD_NECK_x = X(28);
F_HEAD_NECK_y = X(29);
F_HEAD_NECK_z = X(30);
M_Ankle_L_x = X(31);
M_Ankle_L_y = X(32);
M_Ankle_L_z = X(33);
M_Ankle_R_x = X(34);
M_Ankle_R_y = X(35);
M_Ankle_R_z = X(36);
M_Knee_L_x = X(37);
M_Knee_L_y = X(38);
M_Knee_L_z = X(39);
M_Knee_R_x = X(40);
M_Knee_R_y = X(41);
M_Knee_R_z = X(42);
M_Hip_L_x = X(43);
M_Hip_L_y = X(44);
M_Hip_L_z = X(45);
M_Hip_R_x = X(46);
M_Hip_R_y = X(47);
M_Hip_R_z = X(48);
M_Pelvis_x = X(49);
M_Pelvis_y = X(50);
M_Pelvis_z = X(51);
M_Head_Neck_x = X(52);
M_Head_Neck_y = X(53);
M_Head_Neck_z = X(54);
r_1_L_x = X(55);
%r_1_L_y = X(56);
r_1_R_x = X(56);
%r_1_R_y = X(58);




F(1) = F_GRF_L_x - F_ANKLE_L_x - m_foot_L*g_x;
F(2) = F_GRF_L_y - F_ANKLE_L_y - m_foot_L*g_y;
F(3) = F_GRF_L_z - F_ANKLE_L_z - m_foot_L*g_z;

F(4) = F_GRF_R_x - F_ANKLE_R_x - m_foot_R*g_x;
F(5) = F_GRF_R_y - F_ANKLE_R_y - m_foot_R*g_y;
F(6) = F_GRF_R_z - F_ANKLE_R_z - m_foot_R*g_z;

F(7) = F_ANKLE_L_x - F_KNEE_L_x - m_leg_L*g_x;
F(8) = F_ANKLE_L_y - F_KNEE_L_y - m_leg_L*g_y;
F(9) = F_ANKLE_L_z - F_KNEE_L_z - m_leg_L*g_z;

F(10) = F_ANKLE_R_x - F_KNEE_R_x - m_leg_R*g_x;
F(11) = F_ANKLE_R_y - F_KNEE_R_y - m_leg_R*g_y;
F(12) = F_ANKLE_R_z - F_KNEE_R_z - m_leg_R*g_z;

F(13) = F_KNEE_L_x - F_HIP_L_x - m_thigh_L*g_x;
F(14) = F_KNEE_L_y - F_HIP_L_y - m_thigh_L*g_y;
F(15) = F_KNEE_L_z - F_HIP_L_z - m_thigh_L*g_z;

F(16) = F_KNEE_R_x - F_HIP_R_x - m_thigh_R*g_x;
F(17) = F_KNEE_R_y - F_HIP_R_y - m_thigh_R*g_y;
F(18) = F_KNEE_R_z - F_HIP_R_z - m_thigh_R*g_z;

F(19) = F_HIP_L_x + F_HIP_R_x - F_PELVIS_x - m_pelvis*g_x;
F(20) = F_HIP_L_y + F_HIP_R_y - F_PELVIS_y - m_pelvis*g_y;
F(21) = F_HIP_L_z + F_HIP_R_z - F_PELVIS_z - m_pelvis*g_z;

F(22) = F_PELVIS_x - F_HEAD_NECK_x - m_torso_arms*g_x;
F(23) = F_PELVIS_y - F_HEAD_NECK_y - m_torso_arms*g_y;
F(24) = F_PELVIS_z - F_HEAD_NECK_z - m_torso_arms*g_z;

F(25) = r_1_L_y*F_GRF_L_z - r_1_L_z*F_GRF_L_y...
    - M_Ankle_L_x - r_2_L_y*F_ANKLE_L_z + r_2_L_z*F_ANKLE_L_y;
F(26) = r_1_L_z*F_GRF_L_x - r_1_L_x*F_GRF_L_z...
    - M_Ankle_L_y - r_2_L_z*F_ANKLE_L_x + r_2_L_x*F_ANKLE_L_z;
F(27) =  r_1_L_x*F_GRF_L_y - r_1_L_y*F_GRF_L_x...
    - M_Ankle_L_z - r_2_L_x*F_ANKLE_L_y + r_2_L_y*F_ANKLE_L_x;

F(28) = r_1_R_y*F_GRF_R_z - r_1_R_z*F_GRF_R_y...
    - M_Ankle_R_x - r_2_R_y*F_ANKLE_R_z + r_2_R_z*F_ANKLE_R_y;
F(29) = r_1_R_z*F_GRF_R_x - r_1_R_x*F_GRF_R_z...
    - M_Ankle_R_y - r_2_R_z*F_ANKLE_R_x + r_2_R_x*F_ANKLE_R_z;
F(30) =  r_1_R_x*F_GRF_R_y - r_1_R_y*F_GRF_R_x...
    - M_Ankle_R_z - r_2_R_x*F_ANKLE_R_y + r_2_R_y*F_ANKLE_R_x;

F(31) = M_Ankle_L_x - M_Knee_L_x + r_3_L_y*F_ANKLE_L_z ...
    - r_3_L_z*F_ANKLE_L_y - r_4_L_y*F_KNEE_L_z ...
    + r_4_L_z*F_KNEE_L_y;

F(32) = M_Ankle_L_y -M_Knee_L_y ...
    + r_3_L_z*F_ANKLE_L_x - r_3_L_x*F_ANKLE_L_z...
    - r_4_L_z*F_KNEE_L_x + r_4_L_x*F_KNEE_L_z;

F(33) = M_Ankle_L_z -M_Knee_L_z ...
    + r_3_L_x*F_ANKLE_L_y - r_3_L_y*F_ANKLE_L_x...
    - r_4_L_x*F_KNEE_L_y + r_4_L_y*F_KNEE_L_x;

%Equation 7 R
F(34) = M_Ankle_R_x -M_Knee_R_x ...
    + r_3_R_y*F_ANKLE_R_z - r_3_R_z*F_ANKLE_R_y...
    - r_4_R_y*F_KNEE_R_z + r_4_R_z*F_KNEE_R_y;

F(35) = M_Ankle_R_y -M_Knee_R_y ...
    + r_3_R_z*F_ANKLE_R_x - r_3_R_x*F_ANKLE_R_z...
    - r_4_R_z*F_KNEE_R_x + r_4_R_x*F_KNEE_R_z;

F(36) = M_Ankle_R_z -M_Knee_R_z ...
    + r_3_R_x*F_ANKLE_R_y - r_3_R_y*F_ANKLE_R_x...
    - r_4_R_x*F_KNEE_R_y + r_4_R_y*F_KNEE_R_x;

F(37) = M_Knee_L_x - M_Hip_L_x ...
    + r_5_L_y*F_KNEE_L_z - r_5_L_z*F_KNEE_L_y...
    - r_6_L_y*F_HIP_L_z + r_6_L_z*F_HIP_L_y;
F(38) = M_Knee_L_y - M_Hip_L_y ...
    + r_5_L_z*F_KNEE_L_x - r_5_L_x*F_KNEE_L_z...
    - r_6_L_z*F_HIP_L_x + r_6_L_x*F_HIP_L_z;
F(39) = M_Knee_L_z - M_Hip_L_z ...
    + r_5_L_x*F_KNEE_L_y - r_5_L_y*F_KNEE_L_x...
    - r_6_L_x*F_HIP_L_y + r_6_L_y*F_HIP_L_x;

F(40) = M_Knee_R_x - M_Hip_R_x ...
    + r_5_R_y*F_KNEE_R_z - r_5_R_z*F_KNEE_R_y...
    - r_6_R_y*F_HIP_R_z + r_6_R_z*F_HIP_R_y;
F(41) = M_Knee_R_y - M_Hip_R_y ...
    + r_5_R_z*F_KNEE_R_x - r_5_R_x*F_KNEE_R_z...
    - r_6_R_z*F_HIP_R_x + r_6_R_x*F_HIP_R_z;
F(42) = M_Knee_R_z - M_Hip_R_z ...
    + r_5_R_x*F_KNEE_R_y - r_5_R_y*F_KNEE_R_x...
    - r_6_R_x*F_HIP_R_y + r_6_R_y*F_HIP_R_x;


F(43) = M_Hip_L_x + M_Hip_R_x - M_Pelvis_x...
    + r_7_L_y*F_HIP_L_z - r_7_L_z*F_HIP_L_y...
    - r_7_R_y*F_HIP_R_z + r_7_R_z*F_HIP_R_y;

F(44) = M_Hip_L_y + M_Hip_R_y - M_Pelvis_y...
    + r_7_L_z*F_HIP_L_x - r_7_L_x*F_HIP_L_z...
    - r_7_R_z*F_HIP_R_x + r_7_R_x*F_HIP_R_z;

F(45) = M_Hip_L_z + M_Hip_R_z - M_Pelvis_z...
    + r_7_L_x*F_HIP_L_y - r_7_L_y*F_HIP_L_x...
    - r_7_R_x*F_HIP_R_y + r_7_R_y*F_HIP_R_x;


F(46) = M_Pelvis_x - M_Head_Neck_x...
    + r_8_y*F_PELVIS_z - r_8_z*F_PELVIS_y...
    - r_9_y*F_HEAD_NECK_z + r_9_z*F_HEAD_NECK_y;

F(47) = M_Pelvis_y - M_Head_Neck_y...
    + r_8_z*F_PELVIS_x - r_8_x*F_PELVIS_z...
    - r_9_z*F_HEAD_NECK_x + r_9_x*F_HEAD_NECK_z;

F(48) = M_Pelvis_z - M_Head_Neck_z...
     + r_8_x*F_PELVIS_y - r_8_y*F_PELVIS_x...
    - r_9_x*F_HEAD_NECK_y + r_9_y*F_HEAD_NECK_x;

F(49) = F_HEAD_NECK_z - m_head_neck*g_x;
F(50) = F_HEAD_NECK_y - m_head_neck*g_y;
F(51) = F_HEAD_NECK_z - m_head_neck*g_z;

F(52) = M_Head_Neck_x + r_10_y*F_HEAD_NECK_z - r_10_z*F_HEAD_NECK_y;
F(53) = M_Head_Neck_y + r_10_z*F_HEAD_NECK_x - r_10_x*F_HEAD_NECK_z;
F(54) = M_Head_Neck_z + r_10_x*F_HEAD_NECK_y - r_10_y*F_HEAD_NECK_x;

















