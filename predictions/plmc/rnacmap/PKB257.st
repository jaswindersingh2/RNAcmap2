#Name: PKB257
#Length:  66 
#PageNumber: 2
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
GGGAAACUCCCCGGCCCCGCUGUAGGGGGACCUUCAGCGACAGGGCCAGAACGAAUAAGGUCCCCA
......[(..]....[.......(.(.].............))(.............)...)....
EEEEEEESMMMMMMMMMMMMMMMSBSHHHHHHHHHHHHHHHSSSHHHHHHHHHHHHHSMMMSEEEE
NNNNNNKNNNKNNNNKNNNNNNNNNNNKNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
S1 8..8 "U" 62..62 "C"
S2 24..24 "A" 43..43 "G"
S3 26..26 "G" 42..42 "A"
S4 44..44 "G" 58..58 "A"
H1 27..41 "GGGACCUUCAGCGAC" (26,42) G:A PK{2}
H2 45..57 "GCCAGAACGAAUA" (44,58) G:A 
B1 25..25 "G" (24,43) A:G (26,42) G:A 
M1.1 9..23 "CCCCGGCCCCGCUGU" (8,62) U:C (24,43) A:G PK{1,2}
M1.2 44..43 "" (43,24) G:A (44,58) G:A 
M1.3 59..61 "GGU" (58,44) A:G (62,8) C:U 
E1 1..7 "GGGAAAC" PK{1}
E2 63..66 "CCCA" 
PK1 1bp 7..7 11..11 E1 1..7 M1.1 9..23
PK2 1bp 16..16 28..28 M1.1 9..23 H1 27..41
PK1.1 7 C 11 C
PK2.1 16 C 28 G
NCBP1 24 A 43 G S2
NCBP2 8 U 62 C S1
NCBP3 44 G 58 A S4
NCBP4 26 G 42 A S3
NCBP5 7 C 11 C PK1.1
segment1 1bp 8..8 U 62..62 C
segment2 2bp 24..26 AGG 42..43 AG
segment3 1bp 44..44 G 58..58 A
