#Name: 3eph_E
#Length:  69 
#PageNumber: 3
CUCGUAUGGCGCAGUGGUAGCGCAGCAGAUUGCAAAUCUGUUGGUCCUUAGUUCGAUCCUGAGUGCGAG
(((((((.(((..[.....[))).(((((.......)))))..]{].(((.....}...))))))))))
SSSSSSSMSSSHHHHHHHHHSSSMSSSSSHHHHHHHSSSSSMMMMMMSSSHHHHHHHHHSSSSSSSSSS
NNNNNNNNNNNNNKNNNNNKNNNNNNNNNNNNNNNNNNNNNNNKKKNNNNNNNNNKNNNNNNNNNNNNN
S1 1..7 "CUCGUAU" 63..69 "GUGCGAG"
S2 9..11 "GCG" 21..23 "CGC"
S3 25..29 "GCAGA" 37..41 "UCUGU"
S4 48..50 "UUA" 60..62 "UGA"
H1 12..20 "CAGUGGUAG" (11,21) G:C PK{1,2}
H2 30..36 "UUGCAAA" (29,37) A:U 
H3 51..59 "GUUCGAUCC" (50,60) A:U PK{3}
M1.1 8..8 "G" (7,63) U:G (9,23) G:C 
M1.2 24..24 "A" (23,9) C:G (25,41) G:U 
M1.3 42..47 "UGGUCC" (41,25) U:G (48,62) U:A PK{1,2,3}
M1.4 63..62 "" (62,48) A:U (63,7) G:U 
PK1 1bp 14..14 46..46 H1 12..20 M1.3 42..47
PK2 1bp 20..20 44..44 H1 12..20 M1.3 42..47
PK3 1bp 45..45 56..56 M1.3 42..47 H3 51..59
PK1.1 14 G 46 C
PK2.1 20 G 44 G
PK3.1 45 U 56 A
NCBP1 20 G 44 G PK2.1
segment1 7bp 1..7 CUCGUAU 63..69 GUGCGAG
segment2 3bp 9..11 GCG 21..23 CGC
segment3 5bp 25..29 GCAGA 37..41 UCUGU
segment4 3bp 48..50 UUA 60..62 UGA
