#Name: 4pqv_A
#Length:  68 
#PageNumber: 3
GGGUCAGAUCGGCGAAAGUCGCCACUUCGCCGAGGAGUGCAAUCUGUGAGGCCCCAGGAGGACUGGGU
[{.....((....(.....).....(](..[...))....}))......]..................
EEEEEEESSMMMMSHHHHHSMMMMMSBSHHHHHHSSMMMMMSSEEEEEEEEEEEEEEEEEEEEEEEEE
KKNNNNNNNNNNNNNNNNNNNNNNNNKNNNKNNNNNNNNNKNNNNNNNNKNNNNNNNNNNNNNNNNNN
S1 8..9 "AU" 42..43 "AU"
S2 14..14 "G" 20..20 "C"
S3 26..26 "U" 36..36 "A"
S4 28..28 "C" 35..35 "G"
H1 15..19 "AAAGU" (14,20) G:C 
H2 29..34 "GCCGAG" (28,35) C:G PK{3}
B1 27..27 "U" (26,36) U:A (28,35) C:G PK{1}
M1.1 10..13 "CGGC" (9,42) U:A (14,20) G:C 
M1.2 21..25 "GCCAC" (20,14) C:G (26,36) U:A 
M1.3 37..41 "GUGCA" (36,26) A:U (42,9) A:U PK{2}
E1 1..7 "GGGUCAG" PK{1,2}
E2 44..68 "CUGUGAGGCCCCAGGAGGACUGGGU" PK{3}
PK1 1bp 1..1 27..27 E1 1..7 B1 27..27
PK2 1bp 2..2 41..41 E1 1..7 M1.3 37..41
PK3 1bp 31..31 50..50 H2 29..34 E2 44..68
PK1.1 1 G 27 U
PK2.1 2 G 41 A
PK3.1 31 C 50 G
NCBP1 2 G 41 A PK2.1
segment1 2bp 8..9 AU 42..43 AU
segment2 1bp 14..14 G 20..20 C
segment3 2bp 26..28 UUC 35..36 GA
