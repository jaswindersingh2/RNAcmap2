#Name: PKB130
#Length:  40 
#PageNumber: 3
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
AGGGGUAAAACCCCUCGCCUACGUAAGCGUUAUUACGCCC
.(..[(...{.]))}................(...)....
ESBBBSHHHHHHSSXXXXXXXXXXXXXXXXXSHHHSEEEE
NNNNKNNNNKNKNNKNNNNNNNNNNNNNNNNNNNNNNNNN
S1 2..2 "G" 14..14 "C"
S2 6..6 "U" 13..13 "C"
S3 32..32 "A" 36..36 "C"
H1 7..12 "AAAACC" (6,13) U:C PK{1,2}
H2 33..35 "UUA" (32,36) A:C 
B1 3..5 "GGG" (2,14) G:C (6,13) U:C PK{1}
X1 15..31 "UCGCCUACGUAAGCGUU" (14,2) C:G (32,36) A:C PK{2}
E1 1..1 "A" 
E2 37..40 "GCCC" 
PK1 1bp 5..5 12..12 B1 3..5 H1 7..12
PK2 1bp 10..10 15..15 H1 7..12 X1 15..31
PK1.1 5 G 12 C
PK2.1 10 A 15 U
NCBP1 6 U 13 C S2
NCBP2 32 A 36 C S3
segment1 2bp 2..6 GGGGU 13..14 CC
segment2 1bp 32..32 A 36..36 C
