#Name: PKB330
#Length:  64 
#PageNumber: 2
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
GGGGGCCAUUGAAGGAGCGUUCACGUCGCGGUCCCUGUCAGAUGAAAAUCUGCGAAUCCUUCAA
.[.[.[...(](...[..(......]).(...].]...(((((....))))).)......).).
EEEEEEEEESISMMMMMMSHHHHHHHSMSIIIIIIIIISSSSSHHHHSSSSSISMMMMMMSISE
NKNKNKNNNNKNNNNKNNNNNNNNNKNNNNNNKNKNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
S1 10..10 "U" 63..63 "A"
S2 12..12 "A" 61..61 "U"
S3 19..19 "G" 27..27 "C"
S4 29..29 "C" 54..54 "G"
S5 39..43 "CAGAU" 48..52 "AUCUG"
H1 20..26 "UUCACGU" (19,27) G:C PK{3}
H2 44..47 "GAAA" (43,48) U:A 
I1.1 11..11 "G" (10,63) U:A PK{2}
I1.2 62..62 "C" (61,12) U:A 
I2.1 30..38 "GGUCCCUGU" (29,54) C:G PK{1}
I2.2 53..53 "C" (52,39) G:C 
M1.1 13..18 "AGGAGC" (12,61) A:U (19,27) G:C PK{3}
M1.2 28..28 "G" (27,19) C:G (29,54) C:G 
M1.3 55..60 "AAUCCU" (54,29) G:C (61,12) U:A 
E1 1..9 "GGGGGCCAU" PK{1,2}
E2 64..64 "A" 
PK1 2bp 2..4 33..35 E1 1..9 I2.1 30..38
PK2 1bp 6..6 11..11 E1 1..9 I1.1 11..11
PK3 1bp 16..16 26..26 M1.1 13..18 H1 20..26
PK1.1 2 G 35 C
PK1.2 4 G 33 C
PK2.1 6 C 11 G
PK3.1 16 A 26 U
segment1 2bp 10..12 UGA 61..63 UCA
segment2 1bp 19..19 G 27..27 C
segment3 6bp 29..43 CGGUCCCUGUCAGAU 48..54 AUCUGCG
