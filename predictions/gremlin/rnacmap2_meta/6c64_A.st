#Name: 6c64_A
#Length: ,36,
#PageNumber: 2
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
GCGUACGAAGGAGAGGAGAGGAAGAGGAGAGUACGC
((([.(..[............(...]..)].).)))
SSSIISIIIIIIIIIIIIIIISHHHHHHSIISISSS
NNNKNNNNKNNNNNNNNNNNNNNNNKNNNKNNNNNN
S1 1..3 "GCG" 34..36 "CGC"
S2 6..6 "C" 32..32 "U"
S3 22..22 "A" 29..29 "G"
H1 23..28 "AGAGGA" (22,29) A:G PK{2}
I1.1 4..5 "UA" (3,34) G:C PK{1}
I1.2 33..33 "A" (32,6) U:C 
I2.1 7..21 "GAAGGAGAGGAGAGG" (6,32) C:U PK{2}
I2.2 30..31 "AG" (29,22) G:A PK{1}
PK1 1bp 4..4 30..30 I1.1 4..5 I2.2 30..31
PK2 1bp 9..9 26..26 I2.1 7..21 H1 23..28
PK1.1 4 U 30 A
PK2.1 9 A 26 G
NCBP1 22 A 29 G S3
NCBP2 6 C 32 U S2
NCBP3 9 A 26 G PK2.1
segment1 5bp 1..22 GCGUACGAAGGAGAGGAGAGGA 29..36 GAGUACGC
