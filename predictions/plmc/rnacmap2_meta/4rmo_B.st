#Name: 4rmo_B
#Length:  45 
#PageNumber: 3
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
AAGUUUACCACUGACCGAUAUGUGGUAUAUAAAUGGUCGGGUUGA
[...(]..[(......{...].)..)......}(.((...)).).
EEEESIIIISHHHHHHHHHHHHSIISXXXXXXXSISSHHHSSISE
KNNNNKNNKNNNNNNNKNNNKNNNNNNNNNNNKNNNNNNNNNNNN
S1 5..5 "U" 26..26 "U"
S2 10..10 "A" 23..23 "U"
S3 34..34 "U" 44..44 "G"
S4 36..37 "GU" 41..42 "GU"
H1 11..22 "CUGACCGAUAUG" (10,23) A:U PK{2,3}
H2 38..40 "CGG" (37,41) U:G 
I1.1 6..9 "UACC" (5,26) U:U PK{1,2}
I1.2 24..25 "GG" (23,10) U:A 
I2.1 35..35 "G" (34,44) U:G 
I2.2 43..43 "U" (42,36) U:G 
X1 27..33 "AUAUAAA" (26,5) U:U (34,44) U:G PK{3}
E1 1..4 "AAGU" PK{1}
E2 45..45 "A" 
PK1 1bp 1..1 6..6 E1 1..4 I1.1 6..9
PK2 1bp 9..9 21..21 I1.1 6..9 H1 11..22
PK3 1bp 17..17 33..33 H1 11..22 X1 27..33
PK1.1 1 A 6 U
PK2.1 9 C 21 U
PK3.1 17 G 33 A
NCBP1 5 U 26 U S1
NCBP2 9 C 21 U PK2.1
NCBP3 17 G 33 A PK3.1
segment1 2bp 5..10 UUACCA 23..26 UGGU
segment2 3bp 34..37 UGGU 41..44 GUUG
