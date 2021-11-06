#Name: 6dtd_C
#Length: ,37,
#PageNumber: 3
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
UGUUGCAUCUGCCUUCUUUUUGAAAGGUAAAAACAAC
.(.((..[.({(.(.(.]...)..)}).)....)).)
ESISSIIIISISISISHHHHHSIISISISIIIISSIS
NNNNNNNKNNKNNNNNNKNNNNNNNKNNNNNNNNNNN
S1 2..2 "G" 37..37 "C"
S2 4..5 "UG" 34..35 "CA"
S3 10..10 "U" 29..29 "A"
S4 12..12 "C" 27..27 "G"
S5 14..14 "U" 25..25 "A"
S6 16..16 "C" 22..22 "G"
H1 17..21 "UUUUU" (16,22) C:G PK{1}
I1.1 3..3 "U" (2,37) G:C 
I1.2 36..36 "A" (35,4) A:U 
I2.1 6..9 "CAUC" (5,34) G:C PK{1}
I2.2 30..33 "AAAA" (29,10) A:U 
I3.1 11..11 "G" (10,29) U:A PK{2}
I3.2 28..28 "U" (27,12) G:C 
I4.1 13..13 "C" (12,27) C:G 
I4.2 26..26 "G" (25,14) A:U PK{2}
I5.1 15..15 "U" (14,25) U:A 
I5.2 23..24 "AA" (22,16) G:C 
E1 1..1 "U" 
PK1 1bp 8..8 18..18 I2.1 6..9 H1 17..21
PK2 1bp 11..11 26..26 I3.1 11..11 I4.2 26..26
PK1.1 8 U 18 U
PK2.1 11 G 26 G
NCBP1 8 U 18 U PK1.1
NCBP2 11 G 26 G PK2.1
segment1 7bp 2..16 GUUGCAUCUGCCUUC 22..37 GAAAGGUAAAAACAAC
