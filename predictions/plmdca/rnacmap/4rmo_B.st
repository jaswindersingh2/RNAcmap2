#Name: 4rmo_B
#Length:  45 
#PageNumber: 2
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
AAGUUUACCACUGACCGAUAUGUGGUAUAUAAAUGGUCGGGUUGA
.....[(..(....([...)..)....].(....)...]...)..
EEEEEESMMSIIIISHHHHSIISMMMMMMSHHHHSMMMMMMMSEE
NNNNNKNNNNNNNNNKNNNNNNNNNNNKNNNNNNNNNNKNNNNNN
S1 7..7 "A" 43..43 "U"
S2 10..10 "A" 23..23 "U"
S3 15..15 "C" 20..20 "A"
S4 30..30 "U" 35..35 "G"
H1 16..19 "CGAU" (15,20) C:A PK{2}
H2 31..34 "AAAU" (30,35) U:G 
I1.1 11..14 "CUGA" (10,23) A:U 
I1.2 21..22 "UG" (20,15) A:C 
M1.1 8..9 "CC" (7,43) A:U (10,23) A:U 
M1.2 24..29 "GGUAUA" (23,10) U:A (30,35) U:G PK{2}
M1.3 36..42 "GUCGGGU" (35,30) G:U (43,7) U:A PK{1}
E1 1..6 "AAGUUU" PK{1}
E2 44..45 "GA" 
PK1 1bp 6..6 39..39 E1 1..6 M1.3 36..42
PK2 1bp 16..16 28..28 H1 16..19 M1.2 24..29
PK1.1 6 U 39 G
PK2.1 16 C 28 U
NCBP1 15 C 20 A S3
NCBP2 16 C 28 U PK2.1
segment1 1bp 7..7 A 43..43 U
segment2 2bp 10..15 ACUGAC 20..23 AUGU
segment3 1bp 30..30 U 35..35 G
