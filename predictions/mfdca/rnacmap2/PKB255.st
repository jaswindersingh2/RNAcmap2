#Name: PKB255
#Length:  56 
#PageNumber: 3
CUCUCUAUCAGAAUGGAUGUCUUGCUGUCAUAACAGAUAGAGAAGGUUGUGGCAGA
...[((....{.{.)(].}.[)).((.((((....}...........])))).)).
EEEESSHHHHHHHHSSHHHHHSSXSSISSSSHHHHHHHHHHHHHHHHHSSSSISSE
NNNKNNNNNNKNKNNNKNKNKNNNNNNNNNNNNNNKNNNNNNNNNNNKNNNNNNNN
S1 5..5 "C" 23..23 "U"
S2 6..6 "U" 15..15 "G"
S3 16..16 "G" 22..22 "U"
S4 25..26 "CU" 54..55 "AG"
S5 28..31 "UCAU" 49..52 "GUGG"
H1 7..14 "AUCAGAAU" (6,15) U:G PK{2,3}
H2 17..21 "AUGUC" (16,22) G:U PK{1,3,4}
H3 32..48 "AACAGAUAGAGAAGGUU" (31,49) U:G PK{2,4}
I1.1 27..27 "G" (26,54) U:A 
I1.2 53..53 "C" (52,28) G:U 
M1.1 6..5 "" (5,23) C:U (6,15) U:G 
M1.2 16..15 "" (15,6) G:U (16,22) G:U 
M1.3 23..22 "" (22,16) U:G (23,5) U:C 
X1 24..24 "G" (23,5) U:C (25,55) C:G 
E1 1..4 "CUCU" PK{1}
E2 56..56 "A" 
PK1 1bp 4..4 17..17 E1 1..4 H2 17..21
PK2 1bp 11..11 36..36 H1 7..14 H3 32..48
PK3 1bp 13..13 19..19 H1 7..14 H2 17..21
PK4 1bp 21..21 48..48 H2 17..21 H3 32..48
PK1.1 4 U 17 A
PK2.1 11 G 36 G
PK3.1 13 A 19 G
PK4.1 21 C 48 U
NCBP1 5 C 23 U S1
NCBP2 11 G 36 G PK2.1
NCBP3 13 A 19 G PK3.1
NCBP4 21 C 48 U PK4.1
segment1 1bp 5..5 C 23..23 U
segment2 1bp 6..6 U 15..15 G
segment3 1bp 16..16 G 22..22 U
segment4 6bp 25..31 CUGUCAU 49..55 GUGGCAG
