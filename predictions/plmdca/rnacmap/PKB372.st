#Name: PKB372
#Length:  41 
#PageNumber: 3
CCUCGGGCACGCGCUGAAAGUGCCCGUUGGAAGGAUUGCGA
..[.[((.(.{]......{.)])).......(.}..)}...
EEEEESSISHHHHHHHHHHHSISSXXXXXXXSHHHHSEEEE
NNKNKNNNNNKKNNNNNNKNNKNNNNNNNNNNNKNNNKNNN
S1 6..7 "GG" 23..24 "CC"
S2 9..9 "A" 21..21 "U"
S3 32..32 "A" 37..37 "U"
H1 10..20 "CGCGCUGAAAG" (9,21) A:U PK{2,3,4}
H2 33..36 "GGAU" (32,37) A:U PK{4}
I1.1 8..8 "C" (7,23) G:C 
I1.2 22..22 "G" (21,9) U:A PK{1}
X1 25..31 "CGUUGGA" (24,6) C:G (32,37) A:U 
E1 1..5 "CCUCG" PK{1,2}
E2 38..41 "GCGA" PK{3}
PK1 1bp 3..3 22..22 E1 1..5 I1.2 22..22
PK2 1bp 5..5 12..12 E1 1..5 H1 10..20
PK3 1bp 11..11 38..38 H1 10..20 E2 38..41
PK4 1bp 19..19 34..34 H1 10..20 H2 33..36
PK1.1 3 U 22 G
PK2.1 5 G 12 C
PK3.1 11 G 38 G
PK4.1 19 A 34 G
NCBP1 11 G 38 G PK3.1
NCBP2 19 A 34 G PK4.1
segment1 3bp 6..9 GGCA 21..24 UGCC
segment2 1bp 32..32 A 37..37 U
