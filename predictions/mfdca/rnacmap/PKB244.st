#Name: PKB244
#Length:  55 
#PageNumber: 3
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
UUAAACAAAUUUUCUUAAAAUUUCUGAGGUUUGUUUAUUUCUUUUAUCAGUAAAU
(...(....[{..........)..)...]..(...)....[.(.}].).......
SIIISHHHHHHHHHHHHHHHHSIISMMMMMMSHHHSMMMMMMSHHHHSEEEEEEE
NNNNNNNNNKKNNNNNNNNNNNNNNNNNKNNNNNNNNNNNKNNNKKNNNNNNNNN
S1 1..1 "U" 25..25 "U"
S2 5..5 "A" 22..22 "U"
S3 32..32 "U" 36..36 "U"
S4 43..43 "U" 48..48 "C"
H1 6..21 "CAAAUUUUCUUAAAAU" (5,22) A:U PK{1,2}
H2 33..35 "GUU" (32,36) U:U 
H3 44..47 "UUAU" (43,48) U:C PK{2,3}
I1.1 2..4 "UAA" (1,25) U:U 
I1.2 23..24 "UC" (22,5) U:A 
M1.1 26..31 "GAGGUU" (25,1) U:U (32,36) U:U PK{1}
M1.2 37..42 "AUUUCU" (36,32) U:U (43,48) U:C PK{3}
E1 49..55 "AGUAAAU" 
PK1 1bp 10..10 29..29 H1 6..21 M1.1 26..31
PK2 1bp 11..11 45..45 H1 6..21 H3 44..47
PK3 1bp 41..41 46..46 M1.2 37..42 H3 44..47
PK1.1 10 U 29 G
PK2.1 11 U 45 U
PK3.1 41 C 46 A
NCBP1 32 U 36 U S3
NCBP2 1 U 25 U S1
NCBP3 43 U 48 C S4
NCBP4 11 U 45 U PK2.1
NCBP5 41 C 46 A PK3.1
segment1 2bp 1..5 UUAAA 22..25 UUCU
segment2 1bp 32..32 U 36..36 U
segment3 1bp 43..43 U 48..48 C
