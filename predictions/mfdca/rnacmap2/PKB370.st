#Name: PKB370
#Length:  38 
#PageNumber: 3
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
UGAUGUGCAGCUUCUUAGGUGCAUCUCGUAUAGAAGCG
.[(..(.((...{......))])).[(...]...}..)
EESBBSISSHHHHHHHHHHSSISSXXSHHHHHHHHHHS
NKNNNNNNNNNNKNNNNNNNNKNNNKNNNNKNNNKNNN
S1 3..3 "A" 24..24 "U"
S2 6..6 "U" 23..23 "A"
S3 8..9 "CA" 20..21 "UG"
S4 27..27 "C" 38..38 "G"
H1 10..19 "GCUUCUUAGG" (9,20) A:U PK{2}
H2 28..37 "GUAUAGAAGC" (27,38) C:G PK{2,3}
B1 4..5 "UG" (3,24) A:U (6,23) U:A 
I1.1 7..7 "G" (6,23) U:A 
I1.2 22..22 "C" (21,8) G:C PK{1}
X1 25..26 "CU" (24,3) U:A (27,38) C:G PK{3}
E1 1..2 "UG" PK{1}
PK1 1bp 2..2 22..22 E1 1..2 I1.2 22..22
PK2 1bp 13..13 35..35 H1 10..19 H2 28..37
PK3 1bp 26..26 31..31 X1 25..26 H2 28..37
PK1.1 2 G 22 C
PK2.1 13 U 35 A
PK3.1 26 U 31 U
NCBP1 26 U 31 U PK3.1
segment1 4bp 3..9 AUGUGCA 20..24 UGCAU
segment2 1bp 27..27 C 38..38 G
