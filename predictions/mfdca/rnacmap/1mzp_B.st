#Name: 1mzp_B
#Length: ,55,
#PageNumber: 3
#Warning: Structure contains linked PK-segments of same sizes 2 and 2. Using PK brackets for the more 5' segment
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
GGGAUGCGUAGGAUAGGUGGGAGCCGCAAGGCGCCGGUGAAAUACCACCCUUCCC
....[...[..({...(((....[[((..]].).).....].}..))).]..)..
EEEEEEEEEEESIIIISSSIIIIIISSHHHHHSBSIIIIIIIIIISSSIIIISEE
NNNNKNNNKNNNKNNNNNNNNNNKKNNNNKKNNNNNNNNNKNKNNNNNNKNNNNN
S1 12..12 "G" 53..53 "C"
S2 17..19 "GUG" 46..48 "CAC"
S3 26..26 "G" 35..35 "C"
S4 27..27 "C" 33..33 "G"
H1 28..32 "AAGGC" (27,33) C:G PK{4}
B1 34..34 "C" (33,27) G:C (35,26) C:G 
I1.1 13..16 "AUAG" (12,53) G:C PK{3}
I1.2 49..52 "CCUU" (48,17) C:G PK{1}
I2.1 20..25 "GGAGCC" (19,46) G:C PK{4}
I2.2 36..45 "GGUGAAAUAC" (35,26) C:G PK{2,3}
E1 1..11 "GGGAUGCGUAG" PK{1,2}
E2 54..55 "CC" 
PK1 1bp 5..5 50..50 E1 1..11 I1.2 49..52
PK2 1bp 9..9 41..41 E1 1..11 I2.2 36..45
PK3 1bp 13..13 43..43 I1.1 13..16 I2.2 36..45
PK4 2bp 24..25 30..31 I2.1 20..25 H1 28..32
PK1.1 5 U 50 C
PK2.1 9 U 41 A
PK3.1 13 A 43 U
PK4.1 24 C 31 G
PK4.2 25 C 30 G
NCBP1 5 U 50 C PK1.1
segment1 6bp 12..27 GAUAGGUGGGAGCCGC 33..53 GCCGGUGAAAUACCACCCUUC
