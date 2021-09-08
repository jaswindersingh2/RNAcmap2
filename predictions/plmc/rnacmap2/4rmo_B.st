#Name: 4rmo_B
#Length: ,45,
#PageNumber: 3
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
AAGUUUACCACUGACCGAUAUGUGGUAUAUAAAUGGUCGGGUUGA
......(..[....(.{.....].[....(..].))......).}
EEEEEESIIIIIIISBBBBBBBBBBBBBBSHHHHSSIIIIIISEE
NNNNNNNNNKNNNNNNKNNNNNKNKNNNNNNNKNNNNNNNNNNNK
S1 7..7 "A" 43..43 "U"
S2 15..15 "C" 36..36 "G"
S3 30..30 "U" 35..35 "G"
H1 31..34 "AAAU" (30,35) U:G PK{3}
B1 16..29 "CGAUAUGUGGUAUA" (15,36) C:G (30,35) U:G PK{1,2,3}
I1.1 8..14 "CCACUGA" (7,43) A:U PK{1}
I1.2 37..42 "UCGGGU" (36,15) G:C 
E1 1..6 "AAGUUU" 
E2 44..45 "GA" PK{2}
PK1 1bp 10..10 23..23 I1.1 8..14 B1 16..29
PK2 1bp 17..17 45..45 B1 16..29 E2 44..45
PK3 1bp 25..25 33..33 B1 16..29 H1 31..34
PK1.1 10 A 23 U
PK2.1 17 G 45 A
PK3.1 25 G 33 A
NCBP1 17 G 45 A PK2.1
NCBP2 25 G 33 A PK3.1
segment1 3bp 7..30 ACCACUGACCGAUAUGUGGUAUAU 35..43 GGUCGGGUU
