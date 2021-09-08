#Name: 1g1x_E
#Length:  44 
#PageNumber: 2
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
ACGCCGAUGGCGAAGGCAGCCACCUGGUCCACCCGUGACGCUUU
.[(...(..........)..]..(...)......).........
EESMMMSHHHHHHHHHHSMMMMMSHHHSMMMMMMSEEEEEEEEE
NKNNNNNNNNNNNNNNNNNNKNNNNNNNNNNNNNNNNNNNNNNN
S1 3..3 "G" 35..35 "G"
S2 7..7 "A" 18..18 "A"
S3 24..24 "C" 28..28 "U"
H1 8..17 "UGGCGAAGGC" (7,18) A:A 
H2 25..27 "UGG" (24,28) C:U 
M1.1 4..6 "CCG" (3,35) G:G (7,18) A:A 
M1.2 19..23 "GCCAC" (18,7) A:A (24,28) C:U PK{1}
M1.3 29..34 "CCACCC" (28,24) U:C (35,3) G:G 
E1 1..2 "AC" PK{1}
E2 36..44 "UGACGCUUU" 
PK1 1bp 2..2 21..21 E1 1..2 M1.2 19..23
PK1.1 2 C 21 C
NCBP1 7 A 18 A S2
NCBP2 24 C 28 U S3
NCBP3 3 G 35 G S1
NCBP4 2 C 21 C PK1.1
segment1 1bp 3..3 G 35..35 G
segment2 1bp 7..7 A 18..18 A
segment3 1bp 24..24 C 28..28 U
