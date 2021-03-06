#Name: 4rmo_B
#Length:  45 
#PageNumber: 4
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
AAGUUUACCACUGACCGAUAUGUGGUAUAUAAAUGGUCGGGUUGA
...[{..(.(......<....])..})........(.>...)...
EEEEEEESISHHHHHHHHHHHHSIIISXXXXXXXXSHHHHHSEEE
NNNKKNNNNNNNNNNNKNNNNKNNNKNNNNNNNNNNNKNNNNNNN
S1 8..8 "C" 27..27 "A"
S2 10..10 "A" 23..23 "U"
S3 36..36 "G" 42..42 "U"
H1 11..22 "CUGACCGAUAUG" (10,23) A:U PK{1,3}
H2 37..41 "UCGGG" (36,42) G:U PK{3}
I1.1 9..9 "C" (8,27) C:A 
I1.2 24..26 "GGU" (23,10) U:A PK{2}
X1 28..35 "UAUAAAUG" (27,8) A:C (36,42) G:U 
E1 1..7 "AAGUUUA" PK{1,2}
E2 43..45 "UGA" 
PK1 1bp 4..4 22..22 E1 1..7 H1 11..22
PK2 1bp 5..5 26..26 E1 1..7 I1.2 24..26
PK3 1bp 17..17 38..38 H1 11..22 H2 37..41
PK1.1 4 U 22 G
PK2.1 5 U 26 U
PK3.1 17 G 38 C
NCBP1 8 C 27 A S1
NCBP2 5 U 26 U PK2.1
segment1 2bp 8..10 CCA 23..27 UGGUA
segment2 1bp 36..36 G 42..42 U
