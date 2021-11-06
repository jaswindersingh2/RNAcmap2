#Name: 5mga_B
#Length:  40 
#PageNumber: 3
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
AAUUUCUACUGUUGUAGAUGAGAAGUCAUUUAAUAAGGCC
.[({.])}........(......)................
EESHHHSXXXXXXXXXSHHHHHHSEEEEEEEEEEEEEEEE
NKNKNKNKNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
S1 3..3 "U" 7..7 "U"
S2 17..17 "G" 24..24 "A"
H1 4..6 "UUC" (3,7) U:U PK{1,2}
H2 18..23 "AUGAGA" (17,24) G:A 
X1 8..16 "ACUGUUGUA" (7,3) U:U (17,24) G:A PK{2}
E1 1..2 "AA" PK{1}
E2 25..40 "GUCAUUUAAUAAGGCC" 
PK1 1bp 2..2 6..6 E1 1..2 H1 4..6
PK2 1bp 4..4 8..8 H1 4..6 X1 8..16
PK1.1 2 A 6 C
PK2.1 4 U 8 A
NCBP1 3 U 7 U S1
NCBP2 17 G 24 A S2
NCBP3 2 A 6 C PK1.1
segment1 1bp 3..3 U 7..7 U
segment2 1bp 17..17 G 24..24 A
