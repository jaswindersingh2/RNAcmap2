#Name: 5tpy_A
#Length:  71 
#PageNumber: 3
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
GGGUCAGGCCGGCGAAAGUCGCCACAGUUUGGGGAAAGCUGUGCAGCCUGUAACCCCCCCACGAAAGUGGG
..........((.....).)....((.............)).............[..(..{(..]..}.))
EEEEEEEEEESSHHHHHSBSXXXXSSHHHHHHHHHHHHHSSXXXXXXXXXXXXXXXXSBBBSHHHHHHHSS
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNKNNNNNKNNNKNNKNNN
S1 11..11 "G" 20..20 "C"
S2 12..12 "G" 18..18 "G"
S3 25..26 "CA" 40..41 "UG"
S4 58..58 "C" 71..71 "G"
S5 62..62 "C" 70..70 "G"
H1 13..17 "CGAAA" (12,18) G:G 
H2 27..39 "GUUUGGGGAAAGC" (26,40) A:U 
H3 63..69 "GAAAGUG" (62,70) C:G PK{1,2}
B1 19..19 "U" (18,12) G:G (20,11) C:G 
B2 59..61 "CCA" (58,71) C:G (62,70) C:G PK{2}
X1 21..24 "GCCA" (20,11) C:G (25,41) C:G 
X2 42..57 "UGCAGCCUGUAACCCC" (41,25) G:C (58,71) C:G PK{1}
E1 1..10 "GGGUCAGGCC" 
PK1 1bp 55..55 65..65 X2 42..57 H3 63..69
PK2 1bp 61..61 68..68 B2 59..61 H3 63..69
PK1.1 55 C 65 A
PK2.1 61 A 68 U
NCBP1 12 G 18 G S2
NCBP2 55 C 65 A PK1.1
segment1 2bp 11..12 GG 18..20 GUC
segment2 2bp 25..26 CA 40..41 UG
segment3 2bp 58..62 CCCAC 70..71 GG
