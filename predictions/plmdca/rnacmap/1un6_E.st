#Name: 1un6_E
#Length:  61 
#PageNumber: 2
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
GCCGGCCACACCUACGGGGCCUGGUUAGUACCUGGGAAACCUGGGAAUACCAGGUGCCGGC
([(.....(...)..[...................(..]).................])).
SBSMMMMMSHHHSMMMMMMMMMMMMMMMMMMMMMMSHHHSMMMMMMMMMMMMMMMMMMSSE
NKNNNNNNNNNNNNNKNNNNNNNNNNNNNNNNNNNNNNKNNNNNNNNNNNNNNNNNNKNNN
S1 1..1 "G" 60..60 "G"
S2 3..3 "C" 59..59 "G"
S3 9..9 "C" 13..13 "U"
S4 36..36 "G" 40..40 "C"
H1 10..12 "ACC" (9,13) C:U 
H2 37..39 "AAA" (36,40) G:C PK{2}
B1 2..2 "C" (1,60) G:G (3,59) C:G PK{1}
M1.1 4..8 "GGCCA" (3,59) C:G (9,13) C:U 
M1.2 14..35 "ACGGGGCCUGGUUAGUACCUGG" (13,9) U:C (36,40) G:C PK{2}
M1.3 41..58 "CUGGGAAUACCAGGUGCC" (40,36) C:G (59,3) G:C PK{1}
E1 61..61 "C" 
PK1 1bp 2..2 58..58 B1 2..2 M1.3 41..58
PK2 1bp 16..16 39..39 M1.2 14..35 H2 37..39
PK1.1 2 C 58 C
PK2.1 16 G 39 A
NCBP1 9 C 13 U S3
NCBP2 1 G 60 G S1
NCBP3 2 C 58 C PK1.1
NCBP4 16 G 39 A PK2.1
segment1 2bp 1..3 GCC 59..60 GG
segment2 1bp 9..9 C 13..13 U
segment3 1bp 36..36 G 40..40 C
