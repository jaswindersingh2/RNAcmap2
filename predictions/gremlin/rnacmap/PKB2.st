#Name: PKB2
#Length:  50 
#PageNumber: 4
GGGAAACGGAGUGCGCGGCACCGUCCGCGGAACAAACGGAGAAGGCAGCU
.[......{(<(.......)......].(...........}...>..).)
EEEEEEEEESMSHHHHHHHSMMMMMMMMSHHHHHHHHHHHHHHHHHHSMS
NKNNNNNNKNKNNNNNNNNNNNNNNNKNNNNNNNNNNNNNKNNNKNNNNN
S1 10..10 "A" 50..50 "U"
S2 12..12 "U" 20..20 "A"
S3 29..29 "G" 48..48 "G"
H1 13..19 "GCGCGGC" (12,20) U:A 
H2 30..47 "GAACAAACGGAGAAGGCA" (29,48) G:G PK{2,3}
M1.1 11..11 "G" (10,50) A:U (12,20) U:A PK{3}
M1.2 21..28 "CCGUCCGC" (20,12) A:U (29,48) G:G PK{1}
M1.3 49..49 "C" (48,29) G:G (50,10) U:A 
E1 1..9 "GGGAAACGG" PK{1,2}
PK1 1bp 2..2 27..27 E1 1..9 M1.2 21..28
PK2 1bp 9..9 41..41 E1 1..9 H2 30..47
PK3 1bp 11..11 45..45 M1.1 11..11 H2 30..47
PK1.1 2 G 27 G
PK2.1 9 G 41 G
PK3.1 11 G 45 G
NCBP1 29 G 48 G S3
NCBP2 2 G 27 G PK1.1
NCBP3 9 G 41 G PK2.1
NCBP4 11 G 45 G PK3.1
segment1 1bp 10..10 A 50..50 U
segment2 1bp 12..12 U 20..20 A
segment3 1bp 29..29 G 48..48 G
