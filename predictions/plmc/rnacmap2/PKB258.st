#Name: PKB258
#Length:  60 
#PageNumber: 3
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
GGGAAACGGGAACUGGGCUUGGGACAAGAGCCAUCCCAAGUCCAAGGCCAAGUAGGCUCG
............[...{(.{..(](.}....[.))}..])...............(...)
EEEEEEEEEEEEEEEEESIIIISBSHHHHHHHHSSIIIISXXXXXXXXXXXXXXXSHHHS
NNNNNNNNNNNNKNNNKNNKNNNKNNKNNNNKNNNKNNKNNNNNNNNNNNNNNNNNNNNN
S1 18..18 "C" 40..40 "G"
S2 23..23 "G" 35..35 "C"
S3 25..25 "C" 34..34 "U"
S4 56..56 "G" 60..60 "G"
H1 26..33 "AAGAGCCA" (25,34) C:U PK{3,4}
H2 57..59 "CUC" (56,60) G:G 
B1 24..24 "A" (23,35) G:C (25,34) C:U PK{1}
I1.1 19..22 "UUGG" (18,40) C:G PK{3}
I1.2 36..39 "CCAA" (35,23) C:G PK{2,4}
X1 41..55 "UCCAAGGCCAAGUAG" (40,18) G:C (56,60) G:G 
E1 1..17 "GGGAAACGGGAACUGGG" PK{1,2}
PK1 1bp 13..13 24..24 E1 1..17 B1 24..24
PK2 1bp 17..17 36..36 E1 1..17 I1.2 36..39
PK3 1bp 20..20 27..27 I1.1 19..22 H1 26..33
PK4 1bp 32..32 39..39 H1 26..33 I1.2 36..39
PK1.1 13 C 24 A
PK2.1 17 G 36 C
PK3.1 20 U 27 A
PK4.1 32 C 39 A
NCBP1 56 G 60 G S4
NCBP2 25 C 34 U S3
NCBP3 13 C 24 A PK1.1
NCBP4 32 C 39 A PK4.1
segment1 3bp 18..25 CUUGGGAC 34..40 UCCCAAG
segment2 1bp 56..56 G 60..60 G
