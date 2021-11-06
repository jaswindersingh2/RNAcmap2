#Name: 5o6u_A
#Length:  43 
#PageNumber: 3
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
CUUAGAAAGCAAACCUAUAACCGUUCACCGCCGCACAGGCGGC
[(..])[.(..{](...}))........(...........)..
ESHHHSXXSBBBBSHHHHSSXXXXXXXXSHHHHHHHHHHHSEE
KNNNKNKNNNNKKNNNNKNNNNNNNNNNNNNNNNNNNNNNNNN
S1 2..2 "U" 6..6 "A"
S2 9..9 "G" 20..20 "A"
S3 14..14 "C" 19..19 "A"
S4 29..29 "C" 41..41 "G"
H1 3..5 "UAG" (2,6) U:A PK{1}
H2 15..18 "CUAU" (14,19) C:A PK{3}
H3 30..40 "GCCGCACAGGC" (29,41) C:G 
B1 10..13 "CAAA" (9,20) G:A (14,19) C:A PK{2,3}
X1 7..8 "AA" (6,2) A:U (9,20) G:A PK{2}
X2 21..28 "CCGUUCAC" (20,9) A:G (29,41) C:G 
E1 1..1 "C" PK{1}
E2 42..43 "GC" 
PK1 1bp 1..1 5..5 E1 1..1 H1 3..5
PK2 1bp 7..7 13..13 X1 7..8 B1 10..13
PK3 1bp 12..12 18..18 B1 10..13 H2 15..18
PK1.1 1 C 5 G
PK2.1 7 A 13 A
PK3.1 12 A 18 U
NCBP1 9 G 20 A S2
NCBP2 14 C 19 A S3
NCBP3 7 A 13 A PK2.1
segment1 1bp 2..2 U 6..6 A
segment2 2bp 9..14 GCAAAC 19..20 AA
segment3 1bp 29..29 C 41..41 G
