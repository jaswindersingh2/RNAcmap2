#Name: 6aay_B
#Length:  59 
#PageNumber: 2
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
AAAAAGGGUUUAAAAAAUGAAAGUUGGAACUGCUCUCAUUUUGGAGGGUAAUCACAACA
...[.(............].......(([....)..)....].[(.....].).....)
EEEEESMMMMMMMMMMMMMMMMMMMMSSHHHHHSBBSMMMMMMMSHHHHHHHSMMMMMS
NNNKNNNNNNNNNNNNNNKNNNNNNNNNKNNNNNNNNNNNNKNKNNNNNNKNNNNNNNN
S1 6..6 "G" 59..59 "A"
S2 27..27 "G" 37..37 "C"
S3 28..28 "A" 34..34 "U"
S4 45..45 "A" 53..53 "C"
H1 29..33 "ACUGC" (28,34) A:U PK{2}
H2 46..52 "GGGUAAU" (45,53) A:C PK{3}
B1 35..36 "CU" (34,28) U:A (37,27) C:G 
M1.1 7..26 "GGUUUAAAAAAUGAAAGUUG" (6,59) G:A (27,37) G:C PK{1}
M1.2 38..44 "AUUUUGG" (37,27) C:G (45,53) A:C PK{2,3}
M1.3 54..58 "ACAAC" (53,45) C:A (59,6) A:G 
E1 1..5 "AAAAA" PK{1}
PK1 1bp 4..4 19..19 E1 1..5 M1.1 7..26
PK2 1bp 29..29 42..42 H1 29..33 M1.2 38..44
PK3 1bp 44..44 51..51 M1.2 38..44 H2 46..52
PK1.1 4 A 19 G
PK2.1 29 A 42 U
PK3.1 44 G 51 A
NCBP1 45 A 53 C S4
NCBP2 6 G 59 A S1
NCBP3 4 A 19 G PK1.1
NCBP4 44 G 51 A PK3.1
segment1 1bp 6..6 G 59..59 A
segment2 2bp 27..28 GA 34..37 UCUC
segment3 1bp 45..45 A 53..53 C
